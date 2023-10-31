# Copyright (c) 2008-2023 the MRtrix3 contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Covered Software is provided under this License on an "as is"
# basis, without warranty of any kind, either expressed, implied, or
# statutory, including, without limitation, warranties that the
# Covered Software is free of defects, merchantable, fit for a
# particular purpose or non-infringing.
# See the Mozilla Public License v. 2.0 for more details.
#
# For more details, see http://www.mrtrix.org/.



import glob, itertools, json, math, os, shutil, sys, shlex



def usage(cmdline): #pylint: disable=unused-variable
  from mrtrix3 import app, _version #pylint: disable=no-name-in-module, import-outside-toplevel
  cmdline.set_author('Robert E. Smith (robert.smith@florey.edu.au)')
  cmdline.set_synopsis('Perform diffusion image pre-processing using FSL\'s eddy tool; including inhomogeneity distortion correction using FSL\'s topup tool if possible')
  cmdline.add_description('This script is intended to provide convenience of use of the FSL software tools topup and eddy for performing DWI pre-processing, by encapsulating some of the surrounding image data and metadata processing steps. It is intended to simply these processing steps for most commonly-used DWI acquisition strategies, whilst also providing support for some more exotic acquisitions. The "example usage" section demonstrates the ways in which the script can be used based on the (compulsory) -rpe_* command-line options.')
  cmdline.add_description('More information on use of the dwifslpreproc command can be found at the following link: \nhttps://mrtrix.readthedocs.io/en/' + _version.__tag__ + '/dwi_preprocessing/dwifslpreproc.html')
  cmdline.add_description('Note that the MRtrix3 command dwi2mask will automatically be called to derive a processing mask for the FSL command "eddy", which determines which voxels contribute to the estimation of geometric distortion parameters and possibly also the classification of outlier slices. If FSL command "topup" is used to estimate a susceptibility field, then dwi2mask will be executed on the resuts of running FSL command "applytopup" to the input DWIs; otherwise it will be executed directly on the input DWIs. Alternatively, the -eddy_mask option can be specified in order to manually provide such a processing mask. More information on mask derivation from DWI data can be found at: https://mrtrix.readthedocs.io/en/' + _version.__tag__ + '/dwi_preprocessing/masking.html')
  cmdline.add_description('The "-topup_options" and "-eddy_options" command-line options allow the user to pass desired command-line options directly to the FSL commands topup and eddy. The available options for those commands may vary between versions of FSL; users can interrogate such by querying the help pages of the installed software, and/or the FSL online documentation: (topup) https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide ; (eddy) https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide')
  cmdline.add_description('The script will attempt to run the CUDA version of eddy; if this does not succeed for any reason, or is not present on the system, the CPU version will be attempted instead. By default, the CUDA eddy binary found that indicates compilation against the most recent version of CUDA will be attempted; this can be over-ridden by providing a soft-link "eddy_cuda" within your path that links to the binary you wish to be executed.')
  cmdline.add_description('Note that this script does not perform any explicit registration between images provided to topup via the -se_epi option, and the DWI volumes provided to eddy. In some instances (motion between acquisitions) this can result in erroneous application of the inhomogeneity field during distortion correction. Use of the -align_seepi option is advocated in this scenario, which ensures that the first volume in the series provided to topup is also the first volume in the series provided to eddy, guaranteeing alignment. But a prerequisite for this approach is that the image contrast within the images provided to the -se_epi option must match the b=0 volumes present within the input DWI series: this means equivalent TE, TR and flip angle (note that differences in multi-band factors between two acquisitions may lead to differences in TR).')
  cmdline.add_example_usage('A basic DWI acquisition, where all image volumes are acquired in a single protocol with fixed phase encoding',
                            'dwifslpreproc DWI_in.mif DWI_out.mif -rpe_none -pe_dir ap -readout_time 0.55',
                            'Due to use of a single fixed phase encoding, no EPI distortion correction can be applied in this case.')
  cmdline.add_example_usage('DWIs all acquired with a single fixed phase encoding; but additionally a pair of b=0 images with reversed phase encoding to estimate the inhomogeneity field',
                            'mrcat b0_ap.mif b0_pa.mif b0_pair.mif -axis 3; dwifslpreproc DWI_in.mif DWI_out.mif -rpe_pair -se_epi b0_pair.mif -pe_dir ap -readout_time 0.72 -align_seepi',
                            'Here the two individual b=0 volumes are concatenated into a single 4D image series, and this is provided to the script via the -se_epi option. Note that with the -rpe_pair option used here, which indicates that the SE-EPI image series contains one or more pairs of b=0 images with reversed phase encoding, the FIRST HALF of the volumes in the SE-EPI series must possess the same phase encoding as the input DWI series, while the second half are assumed to contain the opposite phase encoding direction but identical total readout time. Use of the -align_seepi option is advocated as long as its use is valid (more information in the Description section).')
  cmdline.add_example_usage('All DWI directions & b-values are acquired twice, with the phase encoding direction of the second acquisition protocol being reversed with respect to the first',
                            'mrcat DWI_lr.mif DWI_rl.mif DWI_all.mif -axis 3; dwifslpreproc DWI_all.mif DWI_out.mif -rpe_all -pe_dir lr -readout_time 0.66',
                            'Here the two acquisition protocols are concatenated into a single DWI series containing all acquired volumes. The direction indicated via the -pe_dir option should be the direction of phase encoding used in acquisition of the FIRST HALF of volumes in the input DWI series; ie. the first of the two files that was provided to the mrcat command. In this usage scenario, the output DWI series will contain the same number of image volumes as ONE of the acquired DWI series (ie. half of the number in the concatenated series); this is because the script will identify pairs of volumes that possess the same diffusion sensitisation but reversed phase encoding, and perform explicit recombination of those volume pairs in such a way that image contrast in regions of inhomogeneity is determined from the stretched rather than the compressed image.')
  cmdline.add_example_usage('Any acquisition scheme that does not fall into one of the example usages above',
                            'mrcat DWI_*.mif DWI_all.mif -axis 3; mrcat b0_*.mif b0_all.mif -axis 3; dwifslpreproc DWI_all.mif DWI_out.mif -rpe_header -se_epi b0_all.mif -align_seepi',
                            'With this usage, the relevant phase encoding information is determined entirely based on the contents of the relevant image headers, and dwifslpreproc prepares all metadata for the executed FSL commands accordingly. This can therefore be used if the particular DWI acquisition strategy used does not correspond to one of the simple examples as described in the prior examples. This usage is predicated on the headers of the input files containing appropriately-named key-value fields such that MRtrix3 tools identify them as such. In some cases, conversion from DICOM using MRtrix3 commands will automatically extract and embed this information; however this is not true for all scanner vendors and/or software versions. In the latter case it may be possible to manually provide these metadata; either using the -json_import command-line option of dwifslpreproc, or the -json_import or one of the -import_pe_* command-line options of MRtrix3\'s mrconvert command (and saving in .mif format) prior to running dwifslpreproc.')
  cmdline.add_citation('Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2015, 125, 1063-1078', is_external=True)
  cmdline.add_citation('Smith, S. M.; Jenkinson, M.; Woolrich, M. W.; Beckmann, C. F.; Behrens, T. E.; Johansen-Berg, H.; Bannister, P. R.; De Luca, M.; Drobnjak, I.; Flitney, D. E.; Niazy, R. K.; Saunders, J.; Vickers, J.; Zhang, Y.; De Stefano, N.; Brady, J. M. & Matthews, P. M. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 2004, 23, S208-S219', is_external=True)
  cmdline.add_citation('Skare, S. & Bammer, R. Jacobian weighting of distortion corrected EPI data. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 5063', condition='If performing recombination of diffusion-weighted volume pairs with opposing phase encoding directions', is_external=True)
  cmdline.add_citation('Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888', condition='If performing EPI susceptibility distortion correction', is_external=True)
  cmdline.add_citation('Andersson, J. L. R.; Graham, M. S.; Zsoldos, E. & Sotiropoulos, S. N. Incorporating outlier detection and replacement into a non-parametric framework for movement and distortion correction of diffusion MR images. NeuroImage, 2016, 141, 556-572', condition='If including "--repol" in -eddy_options input', is_external=True)
  cmdline.add_citation('Andersson, J. L. R.; Graham, M. S.; Drobnjak, I.; Zhang, H.; Filippini, N. & Bastiani, M. Towards a comprehensive framework for movement and distortion correction of diffusion MR images: Within volume movement. NeuroImage, 2017, 152, 450-466', condition='If including "--mporder" in -eddy_options input', is_external=True)
  cmdline.add_citation('Bastiani, M.; Cottaar, M.; Fitzgibbon, S.P.; Suri, S.; Alfaro-Almagro, F.; Sotiropoulos, S.N.; Jbabdi, S.; Andersson, J.L.R. Automated quality control for within and between studies diffusion MRI data using a non-parametric framework for movement and distortion correction. NeuroImage, 2019, 184, 801-812', condition='If using -eddyqc_text or -eddyqc_all option and eddy_quad is installed', is_external=True)
  cmdline.add_argument('input',  help='The input DWI series to be corrected')
  cmdline.add_argument('output', help='The output corrected image series')
  cmdline.add_argument('-json_import', metavar=('file'), help='Import image header information from an associated JSON file (may be necessary to determine phase encoding information)')
  pe_options = cmdline.add_argument_group('Options for manually specifying the phase encoding of the input DWIs')
  pe_options.add_argument('-pe_dir', metavar=('PE'), help='Manually specify the phase encoding direction of the input series; can be a signed axis number (e.g. -0, 1, +2), an axis designator (e.g. RL, PA, IS), or NIfTI axis codes (e.g. i-, j, k)')
  pe_options.add_argument('-readout_time', metavar=('time'), type=float, help='Manually specify the total readout time of the input series (in seconds)')
  distcorr_options = cmdline.add_argument_group('Options for achieving correction of susceptibility distortions')
  distcorr_options.add_argument('-se_epi', metavar=('image'), help='Provide an additional image series consisting of spin-echo EPI images, which is to be used exclusively by topup for estimating the inhomogeneity field (i.e. it will not form part of the output image series)')
  distcorr_options.add_argument('-align_seepi', action='store_true', help='Achieve alignment between the SE-EPI images used for inhomogeneity field estimation, and the DWIs (more information in Description section)')
  distcorr_options.add_argument('-topup_options', metavar=('" TopupOptions"'), help='Manually provide additional command-line options to the topup command (provide a string within quotation marks that contains at least one space, even if only passing a single command-line option to topup)')
  distcorr_options.add_argument('-topup_files', metavar=('prefix'), help='Provide files generated by prior execution of the FSL "topup" command to be utilised by eddy')
  cmdline.flag_mutually_exclusive_options( [ 'topup_files', 'se_epi' ], False )
  cmdline.flag_mutually_exclusive_options( [ 'topup_files', 'align_seepi' ], False )
  cmdline.flag_mutually_exclusive_options( [ 'topup_files', 'topup_options' ], False )
  eddy_options = cmdline.add_argument_group('Options for affecting the operation of the FSL "eddy" command')
  eddy_options.add_argument('-eddy_mask', metavar=('image'), help='Provide a processing mask to use for eddy, instead of having dwifslpreproc generate one internally using dwi2mask')
  eddy_options.add_argument('-eddy_slspec', metavar=('file'), help='Provide a file containing slice groupings for eddy\'s slice-to-volume registration')
  eddy_options.add_argument('-eddy_options', metavar=('" EddyOptions"'), help='Manually provide additional command-line options to the eddy command (provide a string within quotation marks that contains at least one space, even if only passing a single command-line option to eddy)')
  eddyqc_options = cmdline.add_argument_group('Options for utilising EddyQC')
  eddyqc_options.add_argument('-eddyqc_text', metavar=('directory'), help='Copy the various text-based statistical outputs generated by eddy, and the output of eddy_qc (if installed), into an output directory')
  eddyqc_options.add_argument('-eddyqc_all', metavar=('directory'), help='Copy ALL outputs generated by eddy (including images), and the output of eddy_qc (if installed), into an output directory')
  cmdline.flag_mutually_exclusive_options( [ 'eddyqc_text', 'eddyqc_all' ], False )
  app.add_dwgrad_export_options(cmdline)
  app.add_dwgrad_import_options(cmdline)
  rpe_options = cmdline.add_argument_group('Options for specifying the acquisition phase-encoding design; note that one of the -rpe_* options MUST be provided')
  rpe_options.add_argument('-rpe_none', action='store_true', help='Specify that no reversed phase-encoding image data is being provided; eddy will perform eddy current and motion correction only')
  rpe_options.add_argument('-rpe_pair', action='store_true', help='Specify that a set of images (typically b=0 volumes) will be provided for use in inhomogeneity field estimation only (using the -se_epi option)')
  rpe_options.add_argument('-rpe_all', action='store_true', help='Specify that ALL DWIs have been acquired with opposing phase-encoding')
  rpe_options.add_argument('-rpe_header', action='store_true', help='Specify that the phase-encoding information can be found in the image header(s), and that this is the information that the script should use')
  cmdline.flag_mutually_exclusive_options( [ 'rpe_none', 'rpe_pair', 'rpe_all', 'rpe_header' ], True )
  cmdline.flag_mutually_exclusive_options( [ 'rpe_none', 'se_epi' ], False ) # May still technically provide -se_epi even with -rpe_all
  cmdline.flag_mutually_exclusive_options( [ 'rpe_pair', 'topup_files'] ) # Would involve two separate sources of inhomogeneity field information
  cmdline.flag_mutually_exclusive_options( [ 'se_epi', 'topup_files'] ) # Would involve two separate sources of inhomogeneity field information
  cmdline.flag_mutually_exclusive_options( [ 'rpe_header', 'pe_dir' ], False ) # Can't manually provide phase-encoding direction if expecting it to be in the header
  cmdline.flag_mutually_exclusive_options( [ 'rpe_header', 'readout_time' ], False ) # Can't manually provide readout time if expecting it to be in the header






def execute(): #pylint: disable=unused-variable
  from mrtrix3 import CONFIG, MRtrixError #pylint: disable=no-name-in-module, import-outside-toplevel
  from mrtrix3 import app, fsl, image, matrix, path, phaseencoding, run, utils #pylint: disable=no-name-in-module, import-outside-toplevel

  if utils.is_windows():
    raise MRtrixError('Script cannot run on Windows due to FSL dependency')

  image.check_3d_nonunity(path.from_user(app.ARGS.input, False))

  pe_design = ''
  if app.ARGS.rpe_none:
    pe_design = 'None'
  elif app.ARGS.rpe_pair:
    pe_design = 'Pair'
    if not app.ARGS.se_epi:
      raise MRtrixError('If using the -rpe_pair option, the -se_epi option must be used to provide the spin-echo EPI data to be used by topup')
  elif app.ARGS.rpe_all:
    pe_design = 'All'
  elif app.ARGS.rpe_header:
    pe_design = 'Header'
  else:
    raise MRtrixError('Must explicitly specify phase-encoding acquisition design (even if none)')

  if app.ARGS.align_seepi and not app.ARGS.se_epi:
    raise MRtrixError('-align_seepi option is only applicable when the -se_epi option is also used')

  fsl_path = os.environ.get('FSLDIR', '')
  if not fsl_path:
    raise MRtrixError('Environment variable FSLDIR is not set; please run appropriate FSL configuration script')

  if not pe_design == 'None':
    topup_config_path = os.path.join(fsl_path, 'etc', 'flirtsch', 'b02b0.cnf')
    if not os.path.isfile(topup_config_path):
      raise MRtrixError('Could not find necessary default config file for FSL topup command (expected location: ' + topup_config_path + ')')
    topup_cmd = fsl.exe_name('topup')

  if not fsl.eddy_binary(True) and not fsl.eddy_binary(False):
    raise MRtrixError('Could not find any version of FSL eddy command')
  fsl_suffix = fsl.suffix()
  app.check_output_path(app.ARGS.output)

  # Export the gradient table to the path requested by the user if necessary
  grad_export_option = app.read_dwgrad_export_options()


  eddyqc_path = None
  eddyqc_files = [ 'eddy_parameters', 'eddy_movement_rms', 'eddy_restricted_movement_rms', \
                   'eddy_post_eddy_shell_alignment_parameters', 'eddy_post_eddy_shell_PE_translation_parameters', \
                   'eddy_outlier_report', 'eddy_outlier_map', 'eddy_outlier_n_stdev_map', 'eddy_outlier_n_sqr_stdev_map', \
                   'eddy_movement_over_time' ]
  if app.ARGS.eddyqc_text:
    eddyqc_path = path.from_user(app.ARGS.eddyqc_text, False)
  elif app.ARGS.eddyqc_all:
    eddyqc_path = path.from_user(app.ARGS.eddyqc_all, False)
    eddyqc_files.extend([ 'eddy_outlier_free_data.nii.gz', 'eddy_cnr_maps.nii.gz', 'eddy_residuals.nii.gz' ])
  if eddyqc_path:
    if os.path.exists(eddyqc_path):
      if os.path.isdir(eddyqc_path):
        if any(os.path.exists(os.path.join(eddyqc_path, filename)) for filename in eddyqc_files):
          if app.FORCE_OVERWRITE:
            app.warn('Output eddy QC directory already contains relevant files; these will be overwritten on completion')
          else:
            raise MRtrixError('Output eddy QC directory already contains relevant files (use -force to override)')
      else:
        if app.FORCE_OVERWRITE:
          app.warn('Target for eddy QC output is not a directory; it will be overwritten on completion')
        else:
          raise MRtrixError('Target for eddy QC output exists, and is not a directory (use -force to override)')


  eddy_manual_options = []
  topup_file_userpath = None
  if app.ARGS.eddy_options:
    # Initially process as a list; we'll convert back to a string later
    eddy_manual_options = app.ARGS.eddy_options.strip().split()
    # Check for erroneous usages before we perform any data importing
    if any(entry.startswith('--mask=') for entry in eddy_manual_options):
      raise MRtrixError('Cannot provide eddy processing mask via -eddy_options "--mask=..." as manipulations are required; use -eddy_mask option instead')
    if any(entry.startswith('--slspec=') for entry in eddy_manual_options):
      raise MRtrixError('Cannot provide eddy slice specification file via -eddy_options "--slspec=..." as manipulations are required; use -eddy_slspec option instead')
    if '--resamp=lsr' in eddy_manual_options:
      raise MRtrixError('dwifslpreproc does not currently support least-squares reconstruction; this cannot be simply passed via -eddy_options')
    eddy_topup_entry = [entry for entry in eddy_manual_options if entry.startswith('--topup=')]
    if len(eddy_topup_entry) > 1:
      raise MRtrixError('Input to -eddy_options contains multiple "--topup=" entries')
    if eddy_topup_entry:
      # -topup_files and -se_epi are mutually exclusive, but need to check in case
      #   pre-calculated topup output files were provided this way instead
      if app.ARGS.se_epi:
        raise MRtrixError('Cannot use both -eddy_options "--topup=" and -se_epi')
      topup_file_userpath = path.from_user(eddy_topup_entry[0][len('--topup='):], False)
      eddy_manual_options = [entry for entry in eddy_manual_options if not entry.startswith('--topup=')]


  # Don't import slspec file directly; just make sure it exists
  if app.ARGS.eddy_slspec and not os.path.isfile(path.from_user(app.ARGS.eddy_slspec, False)):
    raise MRtrixError('Unable to find file \"' + app.ARGS.eddy_slspec + '\" provided via -eddy_slspec option')


  # Attempt to find pre-generated topup files before constructing the scratch directory
  topup_input_movpar = None
  topup_input_fieldcoef = None
  if app.ARGS.topup_files:
    if topup_file_userpath:
      raise MRtrixError('Cannot use -topup_files option and also specify "... --topup=<prefix> ..." within content of -eddy_options')
    topup_file_userpath = path.from_user(app.ARGS.topup_files, False)

  execute_applytopup = pe_design != 'None' or topup_file_userpath
  if execute_applytopup:
    applytopup_cmd = fsl.exe_name('applytopup')

  if topup_file_userpath:
    # Find files based on what the user may or may not have specified:
    # - Path to the movement parameters text file
    # - Path to the field coefficients image
    # - Path prefix including the underscore
    # - Path prefix omitting the underscore

    def check_movpar():
      if not os.path.isfile(topup_input_movpar):
        raise MRtrixError('No topup movement parameter file found based on path "' + topup_file_userpath + '" (expected location: ' + topup_input_movpar + ')')

    def find_fieldcoef(fieldcoef_prefix):
      fieldcoef_candidates = glob.glob(fieldcoef_prefix + '_fieldcoef.nii*')
      if not fieldcoef_candidates:
        raise MRtrixError('No topup field coefficient image found based on path "' + topup_file_userpath + '"')
      if len(fieldcoef_candidates) > 1:
        raise MRtrixError('Multiple topup field coefficient images found based on path "' + topup_file_userpath + '": ' + str(fieldcoef_candidates))
      return fieldcoef_candidates[0]

    if os.path.isfile(topup_file_userpath):
      if topup_file_userpath.endswith('_movpar.txt'):
        topup_input_movpar = topup_file_userpath
        topup_input_fieldcoef = find_fieldcoef(topup_file_userpath[:-len('_movpar.txt')])
      elif topup_file_userpath.endswith('_fieldcoef.nii') or topup_file_userpath.endswith('_fieldcoef.nii.gz'):
        topup_input_fieldcoef = topup_file_userpath
        topup_input_movpar = topup_file_userpath
        if topup_input_movpar.endswith('.gz'):
          topup_input_movpar = topup_input_movpar[:-len('.gz')]
        topup_input_movpar = topup_input_movpar[:-len('_fieldcoef.nii')] + '_movpar.txt'
        check_movpar()
      else:
        raise MRtrixError('Unrecognised file "' + topup_file_userpath + '" specified as pre-calculated topup susceptibility field')
    else:
      topup_input_movpar = topup_file_userpath
      if topup_input_movpar[-1] == '_':
        topup_input_movpar = topup_input_movpar[:-1]
      topup_input_movpar += '_movpar.txt'
      check_movpar()
      topup_input_fieldcoef = find_fieldcoef(topup_input_movpar[:-len('_movpar.txt')])


  # Convert all input images into MRtrix format and store in scratch directory first
  app.make_scratch_dir()

  grad_import_option = app.read_dwgrad_import_options()
  json_import_option = ''
  if app.ARGS.json_import:
    json_import_option = ' -json_import ' + path.from_user(app.ARGS.json_import)
  json_export_option = ' -json_export ' + path.to_scratch('dwi.json', True)
  run.command('mrconvert ' + path.from_user(app.ARGS.input) + ' ' + path.to_scratch('dwi.mif') + grad_import_option + json_import_option + json_export_option)
  if app.ARGS.se_epi:
    image.check_3d_nonunity(path.from_user(app.ARGS.se_epi, False))
    run.command('mrconvert ' + path.from_user(app.ARGS.se_epi) + ' ' + path.to_scratch('se_epi.mif'))
  if topup_file_userpath:
    run.function(shutil.copyfile, topup_input_movpar, path.to_scratch('field_movpar.txt', False))
    # Can't run field spline coefficients image through mrconvert:
    #   topup encodes voxel sizes within the three NIfTI intent parameters, and
    #   applytopup requires that these be set, but mrconvert will wipe them
    run.function(shutil.copyfile, topup_input_fieldcoef, path.to_scratch('field_fieldcoef.nii' + ('.gz' if topup_input_fieldcoef.endswith('.nii.gz') else ''), False))
  if app.ARGS.eddy_mask:
    run.command('mrconvert ' + path.from_user(app.ARGS.eddy_mask) + ' ' + path.to_scratch('eddy_mask.mif') + ' -datatype bit')

  app.goto_scratch_dir()


  # Get information on the input images, and check their validity
  dwi_header = image.Header('dwi.mif')
  if not len(dwi_header.size()) == 4:
    raise MRtrixError('Input DWI must be a 4D image')
  dwi_num_volumes = dwi_header.size()[3]
  app.debug('Number of DWI volumes: ' + str(dwi_num_volumes))
  dwi_num_slices = dwi_header.size()[2]
  app.debug('Number of DWI slices: ' + str(dwi_num_slices))
  dwi_pe_scheme = phaseencoding.get_scheme(dwi_header)
  if app.ARGS.se_epi:
    se_epi_header = image.Header('se_epi.mif')
    # This doesn't necessarily apply any more: May be able to combine e.g. a P>>A from -se_epi with an A>>P b=0 image from the DWIs
  #  if not len(se_epi_header.size()) == 4:
  #    raise MRtrixError('File provided using -se_epi option must contain more than one image volume')
    se_epi_pe_scheme = phaseencoding.get_scheme(se_epi_header)
  if 'dw_scheme' not in dwi_header.keyval():
    raise MRtrixError('No diffusion gradient table found')
  grad = dwi_header.keyval()['dw_scheme']
  if not len(grad) == dwi_num_volumes:
    raise MRtrixError('Number of lines in gradient table (' + str(len(grad)) + ') does not match input image (' + str(dwi_num_volumes) + ' volumes); check your input data')


  # Deal with slice timing information for eddy slice-to-volume correction
  slice_encoding_axis = 2
  eddy_mporder = any(s.startswith('--mporder') for s in eddy_manual_options)
  if eddy_mporder:
    if 'SliceEncodingDirection' in dwi_header.keyval():
      slice_encoding_direction = dwi_header.keyval()['SliceEncodingDirection']
      app.debug('Slice encoding direction: ' + slice_encoding_direction)
      if not slice_encoding_direction.startswith('k'):
        raise MRtrixError('DWI header indicates that 3rd spatial axis is not the slice axis; this is not yet compatible with --mporder option in eddy, nor supported in dwifslpreproc')
      slice_encoding_direction = image.axis2dir(slice_encoding_direction)
    else:
      app.console('No slice encoding direction information present; assuming third axis corresponds to slices')
      slice_encoding_direction = [0,0,1]
    slice_encoding_axis = [ index for index, value in enumerate(slice_encoding_direction) if value ][0]
    slice_groups = [ ]
    slice_timing = [ ]
    # Since there's a chance that we may need to pad this info, we can't just copy this file
    #   to the scratch directory...
    if app.ARGS.eddy_slspec:
      try:
        slice_groups = matrix.load_numeric(path.from_user(app.ARGS.eddy_slspec, False), dtype=int)
        app.debug('Slice groups: ' + str(slice_groups))
      except ValueError:
        try:
          slice_timing = matrix.load_numeric(path.from_user(app.ARGS.eddy_slspec, False), dtype=float)
          app.debug('Slice timing: ' + str(slice_timing))
          app.warn('\"slspec\" file provided to FSL eddy is supposed to contain slice indices for slice groups; '
                   'contents of file \"' + app.ARGS.eddy_slspec + '\" appears to instead be slice timings; '
                   'these data have been imported and will be converted to the appropriate format')
          if len(slice_timing) != dwi_num_slices:
            raise MRtrixError('Cannot use slice timing information from file \"' + app.ARGS.eddy_slspec + '\" for slice-to-volume correction: ' # pylint: disable=raise-missing-from
                              'number of entries (' + str(len(slice_timing)) + ') does not match number of slices (' + str(dwi_num_slices) + ')')
        except ValueError:
          raise MRtrixError('Error parsing eddy \"slspec\" file \"' + app.ARGS.eddy_slspec + '\" ' # pylint: disable=raise-missing-from
                            '(please see FSL eddy help page, specifically the --slspec option)')
    else:
      if 'SliceTiming' not in dwi_header.keyval():
        raise MRtrixError('Cannot perform slice-to-volume correction in eddy: '
                          '-eddy_slspec option not specified, and no slice timing information present in input DWI header')
      slice_timing = dwi_header.keyval()['SliceTiming']
      app.debug('Initial slice timing contents from header: ' + str(slice_timing))
      if slice_timing in ['invalid', 'variable']:
        raise MRtrixError('Cannot use slice timing information in image header for slice-to-volume correction: '
                          'data flagged as "' + slice_timing + '"')
      # Fudges necessary to maniupulate nature of slice timing data in cases where
      #   bad JSON formatting has led to the data not being simply a list of floats
      #   (whether from MRtrix3 DICOM conversion or from anything else)
      if isinstance(slice_timing, str):
        slice_timing = slice_timing.split()
      if not isinstance(slice_timing, list):
        raise MRtrixError('Cannot use slice timing information in image header for slice-to-volume correction: '
                          'data is not a list')
      if len(slice_timing) == 1:
        slice_timing = slice_timing[0]
        if not isinstance(slice_timing, list):
          raise MRtrixError('Cannot use slice timing information in image header for slice-to-volume correction: '
                            'unexpected data format')
      if isinstance(slice_timing[0], list):
        if not all(len(entry) == 1 for entry in slice_timing):
          raise MRtrixError('Cannot use slice timing information in image header for slice-to-volume correction: '
                            'data do not appear to be 1D')
        slice_timing = [ entry[0] for entry in slice_timing ]
      if not all(isinstance(entry, float) for entry in slice_timing):
        try:
          slice_timing = [ float(entry) for entry in slice_timing ]
        except ValueError as exception:
          raise MRtrixError('Cannot use slice timing information in image header for slice-to-volume correction: '
                            'data are not numeric') from exception
      app.debug('Re-formatted slice timing contents from header: ' + str(slice_timing))
      if len(slice_timing) != dwi_num_slices:
        raise MRtrixError('Cannot use slice timing information in image header for slice-to-volume correction: '
                          'number of entries (' + str(len(slice_timing)) + ') does not match number of slices (' + str(dwi_header.size()[2]) + ')')
  elif app.ARGS.eddy_slspec:
    app.warn('-eddy_slspec option provided, but "--mporder=" not provided via -eddy_options; '
             'slice specification file not imported as it would not be utilised by eddy')


  # Use new features of dirstat to query the quality of the diffusion acquisition scheme
  # Need to know the mean b-value in each shell, and the asymmetry value of each shell
  # But don't bother testing / warning the user if they're already controlling for this
  if not app.ARGS.eddy_options or not any(s.startswith('--slm=') for s in app.ARGS.eddy_options.split()):
    shell_bvalues = [ int(round(float(value))) for value in image.mrinfo('dwi.mif', 'shell_bvalues').split() ]
    shell_asymmetries = [ float(value) for value in run.command('dirstat dwi.mif -output asym').stdout.splitlines() ]
    # dirstat will skip any b=0 shell by default; therefore for correspondence between
    #   shell_bvalues and shell_symmetry, need to remove any b=0 from the former
    if len(shell_bvalues) == len(shell_asymmetries) + 1:
      shell_bvalues = shell_bvalues[1:]
    elif len(shell_bvalues) != len(shell_asymmetries):
      raise MRtrixError('Number of b-values reported by mrinfo (' + str(len(shell_bvalues)) + ') does not match number of outputs provided by dirstat (' + str(len(shell_asymmetries)) + ')')
    for bvalue, asymmetry in zip(shell_bvalues, shell_asymmetries):
      if asymmetry >= 0.1:
        app.warn('sampling of b=' + str(bvalue) + ' shell is ' + ('strongly' if asymmetry >= 0.4 else 'moderately') + \
                 ' asymmetric; distortion correction may benefit from use of: ' + \
                 '-eddy_options " ... --slm=linear ... "')


  # Since we want to access user-defined phase encoding information regardless of whether or not
  #   such information is present in the header, let's grab it here
  manual_pe_dir = None
  if app.ARGS.pe_dir:
    manual_pe_dir = [ float(i) for i in phaseencoding.direction(app.ARGS.pe_dir) ]
  app.debug('Manual PE direction: ' + str(manual_pe_dir))
  manual_trt = None
  if app.ARGS.readout_time:
    manual_trt = float(app.ARGS.readout_time)
  app.debug('Manual readout time: ' + str(manual_trt))


  # Utilise the b-value clustering algorithm in src/dwi/shells.*
  shell_indices = [ [ int(i) for i in entry.split(',') ] for entry in image.mrinfo('dwi.mif', 'shell_indices').split(' ') ]
  shell_bvalues = [ float(f) for f in image.mrinfo('dwi.mif', 'shell_bvalues').split(' ')]
  bzero_threshold = float(CONFIG.get('BZeroThreshold', 10.0))

  # For each volume index, store the index of the shell to which it is attributed
  #   (this will make it much faster to determine whether or not two volumes belong to the same shell)
  vol2shell = [ -1 ] * dwi_num_volumes
  for index, volumes in enumerate(shell_indices):
    for volume in volumes:
      vol2shell[volume] = index
  assert all(index >= 0 for index in vol2shell)


  def grads_match(one, two):
    # Are the two volumes assigned to different b-value shells?
    if vol2shell[one] != vol2shell[two]:
      return False
    # Does this shell correspond to b=0?
    if shell_bvalues[vol2shell[one]] <= bzero_threshold:
      return True
    # Dot product between gradient directions
    # First, need to check for zero-norm vectors:
    # - If both are zero, skip this check
    # - If one is zero and the other is not, volumes don't match
    # - If neither is zero, test the dot product
    if any(grad[one][0:3]):
      if not any(grad[two][0:3]):
        return False
      dot_product = grad[one][0]*grad[two][0] + grad[one][1]*grad[two][1] + grad[one][2]*grad[two][2]
      if abs(dot_product) < 0.999:
        return False
    elif any(grad[two][0:3]):
      return False
    return True


  # Manually generate a phase-encoding table for the input DWI based on user input
  dwi_manual_pe_scheme = None
  se_epi_manual_pe_scheme = None
  auto_trt = 0.1
  dwi_auto_trt_warning = False
  if manual_pe_dir:

    if manual_trt:
      trt = manual_trt
    else:
      trt = auto_trt
      dwi_auto_trt_warning = True

    # Still construct the manual PE scheme even with 'None' or 'Pair':
    #   there may be information in the header that we need to compare against
    if pe_design == 'None':
      line = list(manual_pe_dir)
      line.append(trt)
      dwi_manual_pe_scheme = [ line ] * dwi_num_volumes
      app.debug('Manual DWI PE scheme for \'None\' PE design: ' + str(dwi_manual_pe_scheme))

    # With 'Pair', also need to construct the manual scheme for SE EPIs
    elif pe_design == 'Pair':
      line = list(manual_pe_dir)
      line.append(trt)
      dwi_manual_pe_scheme = [ line ] * dwi_num_volumes
      app.debug('Manual DWI PE scheme for \'Pair\' PE design: ' + str(dwi_manual_pe_scheme))
      if len(se_epi_header.size()) != 4:
        raise MRtrixError('If using -rpe_pair option, image provided using -se_epi must be a 4D image')
      se_epi_num_volumes = se_epi_header.size()[3]
      if se_epi_num_volumes%2:
        raise MRtrixError('If using -rpe_pair option, image provided using -se_epi must contain an even number of volumes')
      # Assume that first half of volumes have same direction as series;
      #   second half have the opposite direction
      se_epi_manual_pe_scheme = [ line ] * int(se_epi_num_volumes/2)
      line = [ (-i if i else 0.0) for i in manual_pe_dir ]
      line.append(trt)
      se_epi_manual_pe_scheme.extend( [ line ] * int(se_epi_num_volumes/2) )
      app.debug('Manual SEEPI PE scheme for \'Pair\' PE design: ' + str(se_epi_manual_pe_scheme))

    # If -rpe_all, need to scan through grad and figure out the pairings
    # This will be required if relying on user-specified phase encode direction
    # It will also be required at the end of the script for the manual recombination
    # Update: The possible permutations of volume-matched acquisition is limited within the
    #   context of the -rpe_all option. In particular, the potential for having more
    #   than one b=0 volume within each half means that it is not possible to permit
    #   arbitrary ordering of those pairs, since b=0 volumes would then be matched
    #   despite having the same phase-encoding direction. Instead, explicitly enforce
    #   that volumes must be matched between the first and second halves of the DWI data.
    elif pe_design == 'All':
      if dwi_num_volumes%2:
        raise MRtrixError('If using -rpe_all option, input image must contain an even number of volumes')
      grads_matched = [ dwi_num_volumes ] * dwi_num_volumes
      grad_pairs = [ ]
      app.debug('Commencing gradient direction matching; ' + str(dwi_num_volumes) + ' volumes')
      for index1 in range(int(dwi_num_volumes/2)):
        if grads_matched[index1] == dwi_num_volumes: # As yet unpaired
          for index2 in range(int(dwi_num_volumes/2), dwi_num_volumes):
            if grads_matched[index2] == dwi_num_volumes: # Also as yet unpaired
              if grads_match(index1, index2):
                grads_matched[index1] = index2
                grads_matched[index2] = index1
                grad_pairs.append([index1, index2])
                app.debug('Matched volume ' + str(index1) + ' with ' + str(index2) + ': ' + str(grad[index1]) + ' ' + str(grad[index2]))
                break
          else:
            raise MRtrixError('Unable to determine matching reversed phase-encode direction volume for DWI volume ' + str(index1))
      if not len(grad_pairs) == dwi_num_volumes/2:
        raise MRtrixError('Unable to determine complete matching DWI volume pairs for reversed phase-encode combination')
      # Construct manual PE scheme here:
      #   Regardless of whether or not there's a scheme in the header, need to have it:
      #   if there's one in the header, want to compare to the manually-generated one
      dwi_manual_pe_scheme = [ ]
      for index in range(0, dwi_num_volumes):
        line = list(manual_pe_dir)
        if index >= int(dwi_num_volumes/2):
          line = [ (-i if i else 0.0) for i in line ]
        line.append(trt)
        dwi_manual_pe_scheme.append(line)
      app.debug('Manual DWI PE scheme for \'All\' PE design: ' + str(dwi_manual_pe_scheme))

  else: # No manual phase encode direction defined

    if not pe_design == 'Header':
      raise MRtrixError('If not using -rpe_header, phase encoding direction must be provided using the -pe_dir option')



  def scheme_dirs_match(one, two):
    for line_one, line_two in zip(one, two):
      if not line_one[0:3] == line_two[0:3]:
        return False
    return True

  def scheme_times_match(one, two):
    for line_one, line_two in zip(one, two):
      if abs(line_one[3] - line_two[3]) > 5e-3:
        return False
    return True



  # Determine whether or not the phase encoding table generated manually should be used
  #   (possibly instead of a table present in the image header)
  overwrite_dwi_pe_scheme = False
  if dwi_pe_scheme:
    if manual_pe_dir:
      # Compare manual specification to that read from the header;
      #   overwrite & give warning to user if they differ
      # Bear in mind that this could even be the case for -rpe_all;
      #   relying on earlier code having successfully generated the 'appropriate'
      #   PE scheme for the input volume based on the diffusion gradient table
      if not scheme_dirs_match(dwi_pe_scheme, dwi_manual_pe_scheme):
        app.warn('User-defined phase-encoding direction design does not match what is stored in DWI image header; proceeding with user specification')
        overwrite_dwi_pe_scheme = True
    if manual_trt:
      # Compare manual specification to that read from the header
      if not scheme_times_match(dwi_pe_scheme, dwi_manual_pe_scheme):
        app.warn('User-defined total readout time does not match what is stored in DWI image header; proceeding with user specification')
        overwrite_dwi_pe_scheme = True
    if overwrite_dwi_pe_scheme:
      dwi_pe_scheme = dwi_manual_pe_scheme # May be used later for triggering volume recombination
    else:
      dwi_manual_pe_scheme = None # To guarantee that these generated data are never used
  else:
    # Nothing in the header; rely entirely on user specification
    if pe_design == 'Header':
      raise MRtrixError('No phase encoding information found in DWI image header')
    if not manual_pe_dir:
      raise MRtrixError('No phase encoding information provided either in header or at command-line')
    if dwi_auto_trt_warning:
      app.console('Total readout time not provided at command-line; assuming sane default of ' + str(auto_trt))
    dwi_pe_scheme = dwi_manual_pe_scheme # May be needed later for triggering volume recombination

  # This may be required by -rpe_all for extracting b=0 volumes while retaining phase-encoding information
  import_dwi_pe_table_option = ''
  if dwi_manual_pe_scheme:
    phaseencoding.save('dwi_manual_pe_scheme.txt', dwi_manual_pe_scheme)
    import_dwi_pe_table_option = ' -import_pe_table dwi_manual_pe_scheme.txt'


  # Find the index of the first DWI volume that is a b=0 volume
  # This needs to occur at the outermost loop as it is pertinent information
  #   not only for the -align_seepi option, but also for when the -se_epi option
  #   is not provided at all, and the input to topup is extracted solely from the DWIs
  dwi_first_bzero_index = 0
  for line in grad:
    if line[3] <= bzero_threshold:
      break
    dwi_first_bzero_index += 1
  app.debug('Index of first b=0 image in DWIs is ' + str(dwi_first_bzero_index))


  # Deal with the phase-encoding of the images to be fed to topup (if applicable)
  execute_topup = (not pe_design == 'None') and not topup_file_userpath
  overwrite_se_epi_pe_scheme = False
  se_epi_path = 'se_epi.mif'
  dwi_permvols_preeddy_option = ''
  dwi_permvols_posteddy_option = ''
  dwi_bzero_added_to_se_epi = False
  if app.ARGS.se_epi:

    # Newest version of eddy requires that topup field be on the same grid as the eddy input DWI
    if not image.match(dwi_header, se_epi_header, up_to_dim=3):
      app.console('DWIs and SE-EPI images used for inhomogeneity field estimation are defined on different image grids; '
                  'the latter will be automatically re-gridded to match the former')
      new_se_epi_path = 'se_epi_regrid.mif'
      run.command('mrtransform ' + se_epi_path + ' - -reorient_fod no -interp sinc -template dwi.mif | mrcalc - 0.0 -max ' + new_se_epi_path)
      app.cleanup(se_epi_path)
      se_epi_path = new_se_epi_path
      se_epi_header = image.Header(se_epi_path)

    # 3 possible sources of PE information: DWI header, topup image header, command-line
    # Any pair of these may conflict, and any one could be absent

    # Have to switch here based on phase-encoding acquisition design
    if pe_design == 'Pair':
      # Criteria:
      #   * If present in own header, ignore DWI header entirely -
      #     - If also provided at command-line, look for conflict & report
      #     - If not provided at command-line, nothing to do
      #   * If _not_ present in own header:
      #     - If provided at command-line, infer appropriately
      #     - If not provided at command-line, but the DWI header has that information, infer appropriately
      if se_epi_pe_scheme:
        if manual_pe_dir:
          if not scheme_dirs_match(se_epi_pe_scheme, se_epi_manual_pe_scheme):
            app.warn('User-defined phase-encoding direction design does not match what is stored in SE EPI image header; proceeding with user specification')
            overwrite_se_epi_pe_scheme = True
        if manual_trt:
          if not scheme_times_match(se_epi_pe_scheme, se_epi_manual_pe_scheme):
            app.warn('User-defined total readout time does not match what is stored in SE EPI image header; proceeding with user specification')
            overwrite_se_epi_pe_scheme = True
        if overwrite_se_epi_pe_scheme:
          se_epi_pe_scheme = se_epi_manual_pe_scheme
        else:
          se_epi_manual_pe_scheme = None # To guarantee that these data are never used
      else:
        overwrite_se_epi_pe_scheme = True
        se_epi_pe_scheme = se_epi_manual_pe_scheme

    elif pe_design == 'All':
      # Criteria:
      #   * If present in own header:
      #     - Nothing to do
      #   * If _not_ present in own header:
      #     - Don't have enough information to proceed
      #     - Is this too harsh? (e.g. Have rules by which it may be inferred from the DWI header / command-line)
      if not se_epi_pe_scheme:
        raise MRtrixError('If explicitly including SE EPI images when using -rpe_all option, they must come with their own associated phase-encoding information in the image header')

    elif pe_design == 'Header':
      # Criteria:
      #   * If present in own header:
      #       Nothing to do (-pe_dir option is mutually exclusive)
      #   * If _not_ present in own header:
      #       Cannot proceed
      if not se_epi_pe_scheme:
        raise MRtrixError('No phase-encoding information present in SE-EPI image header')
      # If there is no phase encoding contrast within the SE-EPI series,
      #   try combining it with the DWI b=0 volumes, see if that produces some contrast
      # However, this should probably only be permitted if the -align_seepi option is defined
      se_epi_pe_scheme_has_contrast = 'pe_scheme' in se_epi_header.keyval()
      if not se_epi_pe_scheme_has_contrast:
        if app.ARGS.align_seepi:
          app.console('No phase-encoding contrast present in SE-EPI images; will examine again after combining with DWI b=0 images')
          new_se_epi_path = os.path.splitext(se_epi_path)[0] + '_dwibzeros.mif'
          # Don't worry about trying to produce a balanced scheme here
          run.command('dwiextract dwi.mif - -bzero | mrcat - ' + se_epi_path + ' ' + new_se_epi_path + ' -axis 3')
          se_epi_header = image.Header(new_se_epi_path)
          se_epi_pe_scheme_has_contrast = 'pe_scheme' in se_epi_header.keyval()
          if se_epi_pe_scheme_has_contrast:
            app.cleanup(se_epi_path)
            se_epi_path = new_se_epi_path
            se_epi_pe_scheme = phaseencoding.get_scheme(se_epi_header)
            dwi_bzero_added_to_se_epi = True
            # Delay testing appropriateness of the concatenation of these images
            #   (i.e. differences in contrast) to later
          else:
            raise MRtrixError('No phase-encoding contrast present in SE-EPI images, even after concatenating with b=0 images due to -align_seepi option; '
                              'cannot perform inhomogeneity field estimation')
        else:
          raise MRtrixError('No phase-encoding contrast present in SE-EPI images; cannot perform inhomogeneity field estimation')

    if app.ARGS.align_seepi:

      for field_name, description in { 'EchoTime': 'echo time',
                                       'RepetitionTime': 'repetition time',
                                       'FlipAngle': 'flip angle' }.items():
        dwi_value = dwi_header.keyval().get(field_name)
        se_epi_value = se_epi_header.keyval().get(field_name)
        if dwi_value and se_epi_value and dwi_value != se_epi_value:
          app.warn('It appears that the spin-echo EPI images used for inhomogeneity field estimation have a different ' + description + ' to the DWIs being corrected. '
                   'This may cause issues in estimation of the field, as the first DWI b=0 volume will be added to the input series to topup '
                   'due to use of the -align_seepi option.')

      # If we are using the -se_epi option, and hence the input images to topup have not come from the DWIs themselves,
      #   we need to insert the first b=0 DWI volume to the start of the topup input image. Otherwise, the field estimated
      #   by topup will not be correctly aligned with the volumes as they are processed by eddy.
      #
      # However, there's also a code path by which we may have already performed this addition.
      # If we have already apliced the b=0 volumes from the DWI input with the SE-EPI image
      #   (due to the absence of phase-encoding contrast in the SE-EPI series), we don't want to
      #   re-attempt such a concatenation; the fact that the DWI b=0 images were inserted ahead of
      #   the SE-EPI images means the alignment issue should be dealt with.

      if dwi_first_bzero_index == len(grad) and not dwi_bzero_added_to_se_epi:

        app.warn('Unable to find b=0 volume in input DWIs to provide alignment between topup and eddy; script will proceed as though the -align_seepi option were not provided')

      # If b=0 volumes from the DWIs have already been added to the SE-EPI image due to an
      #   absence of phase-encoding contrast in the latter, we don't need to perform the following
      elif not dwi_bzero_added_to_se_epi:

        run.command('mrconvert dwi.mif dwi_first_bzero.mif -coord 3 ' + str(dwi_first_bzero_index) + ' -axes 0,1,2')
        dwi_first_bzero_pe = dwi_manual_pe_scheme[dwi_first_bzero_index] if overwrite_dwi_pe_scheme else dwi_pe_scheme[dwi_first_bzero_index]

        se_epi_pe_sum = [ 0, 0, 0 ]
        se_epi_volume_to_remove = len(se_epi_pe_scheme)
        for index, line in enumerate(se_epi_pe_scheme):
          se_epi_pe_sum = [ i + j for i, j in zip(se_epi_pe_sum, line[0:3]) ]
          if se_epi_volume_to_remove == len(se_epi_pe_scheme) and line[0:3] == dwi_first_bzero_pe[0:3]:
            se_epi_volume_to_remove = index
        new_se_epi_path = os.path.splitext(se_epi_path)[0] + '_firstdwibzero.mif'
        if (se_epi_pe_sum == [ 0, 0, 0 ]) and (se_epi_volume_to_remove < len(se_epi_pe_scheme)):
          app.console('Balanced phase-encoding scheme detected in SE-EPI series; volume ' + str(se_epi_volume_to_remove) + ' will be removed and replaced with first b=0 from DWIs')
          run.command('mrconvert ' + se_epi_path + ' - -coord 3 ' + ','.join([str(index) for index in range(len(se_epi_pe_scheme)) if not index == se_epi_volume_to_remove]) + ' | mrcat dwi_first_bzero.mif - ' + new_se_epi_path + ' -axis 3')
          # Also need to update the phase-encoding scheme appropriately if it's being set manually
          #   (if embedded within the image headers, should be updated through the command calls)
          if se_epi_manual_pe_scheme:
            first_line = list(manual_pe_dir)
            first_line.append(trt)
            new_se_epi_manual_pe_scheme = [ ]
            new_se_epi_manual_pe_scheme.append(first_line)
            for index, entry in enumerate(se_epi_manual_pe_scheme):
              if not index == se_epi_volume_to_remove:
                new_se_epi_manual_pe_scheme.append(entry)
            se_epi_manual_pe_scheme = new_se_epi_manual_pe_scheme
        else:
          if se_epi_pe_sum == [ 0, 0, 0 ] and se_epi_volume_to_remove == len(se_epi_pe_scheme):
            app.console('Phase-encoding scheme of -se_epi image is balanced, but could not find appropriate volume with which to substitute first b=0 volume from DWIs; first b=0 DWI volume will be inserted to start of series, resulting in an unbalanced scheme')
          else:
            app.console('Unbalanced phase-encoding scheme detected in series provided via -se_epi option; first DWI b=0 volume will be inserted to start of series')
          run.command('mrcat dwi_first_bzero.mif ' + se_epi_path + ' ' + new_se_epi_path + ' -axis 3')
          # Also need to update the phase-encoding scheme appropriately
          if se_epi_manual_pe_scheme:
            first_line = list(manual_pe_dir)
            first_line.append(trt)
            se_epi_manual_pe_scheme = [ first_line, se_epi_manual_pe_scheme ]

        # Ended branching based on balanced-ness of PE acquisition scheme within SE-EPI volumes
        app.cleanup(se_epi_path)
        app.cleanup('dwi_first_bzero.mif')
        se_epi_path = new_se_epi_path

      # Ended branching based on:
      # - Detection of first b=0 volume in DWIs; or
      # - Prior merge of SE-EPI and DWI b=0 volumes due to no phase-encoding contrast in SE-EPI

    # Completed checking for presence of -se_epi option

  elif not pe_design == 'None' and not topup_file_userpath: # No SE EPI images explicitly provided: In some cases, can extract appropriate b=0 images from DWI

    # If using 'All' or 'Header', and haven't been given any topup images, need to extract the b=0 volumes from the series,
    #   preserving phase-encoding information while doing so
    # Preferably also make sure that there's some phase-encoding contrast in there...
    # With -rpe_all, need to write inferred phase-encoding to file and import before using dwiextract so that the phase-encoding
    #   of the extracted b=0's is propagated to the generated b=0 series
    run.command('mrconvert dwi.mif' + import_dwi_pe_table_option + ' - | dwiextract - ' + se_epi_path + ' -bzero')
    se_epi_header = image.Header(se_epi_path)

    # If there's no contrast remaining in the phase-encoding scheme, it'll be written to
    #   PhaseEncodingDirection and TotalReadoutTime rather than pe_scheme
    # In this scenario, we will be unable to run topup, or volume recombination
    if 'pe_scheme' not in se_epi_header.keyval():
      if pe_design == 'All':
        raise MRtrixError('DWI header indicates no phase encoding contrast between b=0 images; cannot proceed with volume recombination-based pre-processing')
      app.warn('DWI header indicates no phase encoding contrast between b=0 images; proceeding without inhomogeneity field estimation')
      execute_topup = False
      run.function(os.remove, se_epi_path)
      se_epi_path = None
      se_epi_header = None


  # If the first b=0 volume in the DWIs is in fact not the first volume (i.e. index zero), we're going to
  #   manually place it at the start of the DWI volumes when they are input to eddy, so that the
  #   first input volume to topup and the first input volume to eddy are one and the same.
  # Note: If at a later date, the statistical outputs from eddy are considered (e.g. motion, outliers),
  #   then this volume permutation will need to be taken into account
  if not topup_file_userpath:
    if dwi_first_bzero_index == len(grad):
      app.warn("No image volumes were classified as b=0 by MRtrix3; no permutation of order of DWI volumes can occur " + \
              "(do you need to adjust config file entry BZeroThreshold?)")
    elif dwi_first_bzero_index:
      app.console('First b=0 volume in input DWIs is volume index ' + str(dwi_first_bzero_index) + '; '
                  'this will be permuted to be the first volume (index 0) when eddy is run')
      dwi_permvols_preeddy_option = ' -coord 3 ' + \
                                            str(dwi_first_bzero_index) + \
                                            ',0' + \
                                            (':' + str(dwi_first_bzero_index-1) if dwi_first_bzero_index > 1 else '') + \
                                            (',' + str(dwi_first_bzero_index+1) if dwi_first_bzero_index < dwi_num_volumes-1 else '') + \
                                            (':' + str(dwi_num_volumes-1) if dwi_first_bzero_index < dwi_num_volumes-2 else '')
      dwi_permvols_posteddy_option = ' -coord 3 1' + \
                                            (':' + str(dwi_first_bzero_index) if dwi_first_bzero_index > 1 else '') + \
                                            ',0' + \
                                            (',' + str(dwi_first_bzero_index+1) if dwi_first_bzero_index < dwi_num_volumes-1 else '') + \
                                            (':' + str(dwi_num_volumes-1) if dwi_first_bzero_index < dwi_num_volumes-2 else '')
      app.debug('mrconvert options for axis permutation:')
      app.debug('Pre: ' + str(dwi_permvols_preeddy_option))
      app.debug('Post: ' + str(dwi_permvols_posteddy_option))



  # This may be required when setting up the topup call
  se_epi_manual_pe_table_option = ''
  if se_epi_manual_pe_scheme:
    phaseencoding.save('se_epi_manual_pe_scheme.txt', se_epi_manual_pe_scheme)
    se_epi_manual_pe_table_option = ' -import_pe_table se_epi_manual_pe_scheme.txt'


  # Need gradient table if running dwi2mask after applytopup to derive a brain mask for eddy
  run.command('mrinfo dwi.mif -export_grad_mrtrix grad.b')
  dwi2mask_algo = CONFIG['Dwi2maskAlgorithm']

  eddy_in_topup_option = ''
  dwi_post_eddy_crop_option = ''
  slice_padded = False
  dwi_path = 'dwi.mif'
  if execute_topup:

    # topup will crash if its input image has a spatial dimension with a non-even size;
    #   presumably due to a downsampling by a factor of 2 in a multi-resolution scheme
    # The newest eddy also requires the output from topup and the input DWIs to have the same size;
    #   therefore this restriction applies to the DWIs as well
    # Rather than crop in this case (which would result in a cropped output image),
    #   duplicate the last slice on any problematic axis, and then crop that extra
    #   slice at the output step
    # By this point, if the input SE-EPI images and DWIs are not on the same image grid, the
    #   SE-EPI images have already been re-gridded to DWI image space;
    odd_axis_count = 0
    for axis_size in dwi_header.size()[:3]:
      if int(axis_size%2):
        odd_axis_count += 1
    if odd_axis_count:
      app.console(str(odd_axis_count) + ' spatial ' + ('axes of DWIs have' if odd_axis_count > 1 else 'axis of DWIs has') + ' non-even size; '
                  'this will be automatically padded for compatibility with topup, and the extra slice' + ('s' if odd_axis_count > 1 else '') + ' erased afterwards')
      for axis, axis_size in enumerate(dwi_header.size()[:3]):
        if int(axis_size%2):
          new_se_epi_path = os.path.splitext(se_epi_path)[0] + '_pad' + str(axis) + '.mif'
          run.command('mrconvert ' + se_epi_path + ' -coord ' + str(axis) + ' ' + str(axis_size-1) + ' - | mrcat ' + se_epi_path + ' - ' + new_se_epi_path + ' -axis ' + str(axis))
          app.cleanup(se_epi_path)
          se_epi_path = new_se_epi_path
          new_dwi_path = os.path.splitext(dwi_path)[0] + '_pad' + str(axis) + '.mif'
          run.command('mrconvert ' + dwi_path + ' -coord ' + str(axis) + ' ' + str(axis_size-1) + ' -clear dw_scheme - | mrcat ' + dwi_path + ' - ' + new_dwi_path + ' -axis ' + str(axis))
          app.cleanup(dwi_path)
          dwi_path = new_dwi_path
          dwi_post_eddy_crop_option += ' -coord ' + str(axis) + ' 0:' + str(axis_size-1)
          if axis == slice_encoding_axis:
            slice_padded = True
            dwi_num_slices += 1
            # If we are padding the slice axis, and performing slice-to-volume correction,
            #   then we need to perform the corresponding padding to the slice timing
            if eddy_mporder:
              # At this point in the script, this information may be encoded either within
              #   the slice timing vector (as imported from the image header), or as
              #   slice groups (i.e. in the format expected by eddy). How these data are
              #   stored affects how the padding is performed.
              if slice_timing:
                slice_timing.append(slice_timing[-1])
              elif slice_groups:
                # Can't edit in place when looping through the list
                new_slice_groups = [ ]
                for group in slice_groups:
                  if axis_size-1 in group:
                    group.append(axis_size)
                  new_slice_groups.append(group)
                slice_groups = new_slice_groups


    # Do the conversion in preparation for topup
    run.command('mrconvert ' + se_epi_path + ' topup_in.nii' + se_epi_manual_pe_table_option + ' -strides -1,+2,+3,+4 -export_pe_table topup_datain.txt')
    app.cleanup(se_epi_path)

    # Run topup
    topup_manual_options = ''
    if app.ARGS.topup_options:
      topup_manual_options = ' ' + app.ARGS.topup_options.strip()
    topup_output = run.command(topup_cmd + ' --imain=topup_in.nii --datain=topup_datain.txt --out=field --fout=field_map' + fsl_suffix + ' --config=' + topup_config_path + ' --verbose' + topup_manual_options)
    with open('topup_output.txt', 'wb') as topup_output_file:
      topup_output_file.write((topup_output.stdout + '\n' + topup_output.stderr + '\n').encode('utf-8', errors='replace'))
    if app.VERBOSITY > 1:
      app.console('Output of topup command:')
      sys.stderr.write(topup_output.stdout + '\n' + topup_output.stderr + '\n')

  if execute_applytopup:

    # Apply the warp field to the input image series to get an initial corrected volume estimate
    # applytopup can't receive the complete DWI input and correct it as a whole, because the phase-encoding
    #   details may vary between volumes
    if dwi_manual_pe_scheme:
      run.command('mrconvert ' + dwi_path + import_dwi_pe_table_option + ' - | mrinfo - -export_pe_eddy applytopup_config.txt applytopup_indices.txt')
    else:
      run.command('mrinfo ' + dwi_path + ' -export_pe_eddy applytopup_config.txt applytopup_indices.txt')

    # Call applytopup separately for each unique phase-encoding
    # This should be the most compatible option with more complex phase-encoding acquisition designs,
    #   since we don't need to worry about applytopup performing volume recombination
    # Plus, recombination doesn't need to be optimal; we're only using this to derive a brain mask
    applytopup_image_list = [ ]
    index = 1
    applytopup_config = matrix.load_matrix('applytopup_config.txt')
    applytopup_indices = matrix.load_vector('applytopup_indices.txt', dtype=int)
    applytopup_volumegroups = [ [ index for index, value in enumerate(applytopup_indices) if value == group ] for group in range(1, len(applytopup_config)+1) ]
    app.debug('applytopup_config: ' + str(applytopup_config))
    app.debug('applytopup_indices: ' + str(applytopup_indices))
    app.debug('applytopup_volumegroups: ' + str(applytopup_volumegroups))
    for index, group in enumerate(applytopup_volumegroups):
      prefix = os.path.splitext(dwi_path)[0] + '_pe_' + str(index)
      input_path = prefix + '.nii'
      json_path = prefix + '.json'
      temp_path = prefix + '_applytopup.nii'
      output_path = prefix + '_applytopup.mif'
      run.command('mrconvert ' + dwi_path + ' ' + input_path + ' -coord 3 ' + ','.join(str(value) for value in group) + ' -strides -1,+2,+3,+4 -json_export ' + json_path)
      run.command(applytopup_cmd + ' --imain=' + input_path + ' --datain=applytopup_config.txt --inindex=' + str(index+1) + ' --topup=field --out=' + temp_path + ' --method=jac')
      app.cleanup(input_path)
      temp_path = fsl.find_image(temp_path)
      run.command('mrconvert ' + temp_path + ' ' + output_path + ' -json_import ' + json_path)
      app.cleanup(json_path)
      app.cleanup(temp_path)
      applytopup_image_list.append(output_path)
      index += 1

    # Use the initial corrected volumes to derive a brain mask for eddy
    if not app.ARGS.eddy_mask:

      dwi2mask_out_path = 'dwi2mask_out.mif'
      if len(applytopup_image_list) == 1:
        dwi2mask_in_path = applytopup_image_list[0]
      else:
        dwi2mask_in_path = 'dwi2mask_in.mif'
        run.command('mrcat ' + ' '.join(applytopup_image_list) + ' ' + dwi2mask_in_path + ' -axis 3')
      run.command('dwi2mask ' + dwi2mask_algo + ' ' + dwi2mask_in_path + ' ' + dwi2mask_out_path)
      run.command('maskfilter ' + dwi2mask_out_path + ' dilate - | mrconvert - eddy_mask.nii -datatype float32 -strides -1,+2,+3')
      if len(applytopup_image_list) > 1:
        app.cleanup(dwi2mask_in_path)
      app.cleanup(dwi2mask_out_path)

    app.cleanup(applytopup_image_list)

    eddy_in_topup_option = ' --topup=field'

  else:

    # Generate a processing mask for eddy based on the uncorrected input DWIs
    if not app.ARGS.eddy_mask:
      dwi2mask_out_path = 'dwi2mask_out.mif'
      run.command('dwi2mask ' + dwi2mask_algo + ' ' + dwi_path + ' ' + dwi2mask_out_path)
      run.command('maskfilter ' + dwi2mask_out_path + ' dilate - | mrconvert - eddy_mask.nii -datatype float32 -strides -1,+2,+3')
      app.cleanup(dwi2mask_out_path)


  # Use user supplied mask for eddy instead of one derived from the images using dwi2mask
  if app.ARGS.eddy_mask:
    if image.match('eddy_mask.mif', dwi_path, up_to_dim=3):
      run.command('mrconvert eddy_mask.mif eddy_mask.nii -datatype float32 -stride -1,+2,+3')
    else:
      app.warn('User-provided processing mask for eddy does not match DWI voxel grid; resampling')
      run.command('mrtransform eddy_mask.mif - -template ' + dwi_path + ' -interp linear | '
                  + 'mrthreshold - -abs 0.5 - | '
                  + 'mrconvert - eddy_mask.nii -datatype float32 -stride -1,+2,+3')
    app.cleanup('eddy_mask.mif')

  # Generate the text file containing slice timing / grouping information if necessary
  if eddy_mporder:
    if slice_timing:
      # This list contains, for each slice, the timing offset between acquisition of the
      #   first slice in the volume, and acquisition of that slice
      # Eddy however requires a text file where each row contains those slices that were
      #   acquired with a single readout, in ordered rows from first slice (group)
      #   acquired to last slice (group) acquired
      if sum(slice_encoding_direction) < 0:
        slice_timing = reversed(slice_timing)
      slice_groups = [ [ x[0] for x in g ] for _, g in itertools.groupby(sorted(enumerate(slice_timing), key=lambda x:x[1]), key=lambda x:x[1]) ] #pylint: disable=unused-variable
      app.debug('Slice timing: ' + str(slice_timing))
      app.debug('Resulting slice groups: ' + str(slice_groups))
    # Variable slice_groups may have already been defined in the correct format.
    #   In that instance, there's nothing to do other than write it to file;
    #   UNLESS the slice encoding direction is known to be reversed, in which case
    #   we need to reverse the timings. Would think that this would however be
    #   rare, given it requires that the slspec text file be provided manually but
    #   SliceEncodingDirection to be present.
    elif slice_groups and sum(slice_encoding_direction) < 0:
      new_slice_groups = [ ]
      for group in new_slice_groups:
        new_slice_groups.append([ dwi_num_slices-index for index in group ])
      app.debug('Slice groups reversed due to negative slice encoding direction')
      app.debug('Original: ' + str(slice_groups))
      app.debug('New: ' + str(new_slice_groups))
      slice_groups = new_slice_groups

    matrix.save_numeric('slspec.txt', slice_groups, add_to_command_history=False, fmt='%d')
    eddy_manual_options.append('--slspec=slspec.txt')


  # Revert eddy_manual_options from a list back to a single string
  eddy_manual_options = (' ' + ' '.join(eddy_manual_options)) if eddy_manual_options else ''


  # Prepare input data for eddy
  run.command('mrconvert ' + dwi_path + import_dwi_pe_table_option + dwi_permvols_preeddy_option + ' eddy_in.nii -strides -1,+2,+3,+4 -export_grad_fsl bvecs bvals -export_pe_eddy eddy_config.txt eddy_indices.txt')
  app.cleanup(dwi_path)

  # Run eddy
  # If a CUDA version is in PATH, run that first; if it fails, re-try using the non-CUDA version
  eddy_all_options = '--imain=eddy_in.nii --mask=eddy_mask.nii --acqp=eddy_config.txt --index=eddy_indices.txt --bvecs=bvecs --bvals=bvals' + eddy_in_topup_option + eddy_manual_options + ' --out=dwi_post_eddy --verbose'
  eddy_cuda_cmd = fsl.eddy_binary(True)
  eddy_openmp_cmd = fsl.eddy_binary(False)
  if eddy_cuda_cmd:
    # If running CUDA version, but OpenMP version is also available, don't stop the script if the CUDA version fails
    try:
      eddy_output = run.command(eddy_cuda_cmd + ' ' + eddy_all_options)
    except run.MRtrixCmdError as exception_cuda:
      if not eddy_openmp_cmd:
        raise
      with open('eddy_cuda_failure_output.txt', 'wb') as eddy_output_file:
        eddy_output_file.write(str(exception_cuda).encode('utf-8', errors='replace'))
      app.console('CUDA version of \'eddy\' was not successful; attempting OpenMP version')
      try:
        eddy_output = run.command(eddy_openmp_cmd + ' ' + eddy_all_options)
      except run.MRtrixCmdError as exception_openmp:
        with open('eddy_openmp_failure_output.txt', 'wb') as eddy_output_file:
          eddy_output_file.write(str(exception_openmp).encode('utf-8', errors='replace'))
        # Both have failed; want to combine error messages
        eddy_cuda_header = ('=' * len(eddy_cuda_cmd)) \
                           + '\n' \
                           + eddy_cuda_cmd \
                           + '\n' \
                           + ('=' * len(eddy_cuda_cmd)) \
                           + '\n'
        eddy_openmp_header = ('=' * len(eddy_openmp_cmd)) \
                             + '\n' \
                             + eddy_openmp_cmd \
                             + '\n' \
                             + ('=' * len(eddy_openmp_cmd)) \
                             + '\n'
        exception_stdout = eddy_cuda_header \
                           + exception_cuda.stdout \
                           + '\n\n' \
                           + eddy_openmp_header \
                           + exception_openmp.stdout \
                           + '\n\n'
        exception_stderr = eddy_cuda_header \
                           + exception_cuda.stderr \
                           + '\n\n' \
                           + eddy_openmp_header \
                           + exception_openmp.stderr \
                           + '\n\n'
        raise run.MRtrixCmdError('eddy* ' + eddy_all_options,
                                 1,
                                 exception_stdout,
                                 exception_stderr)

  else:
    eddy_output = run.command(eddy_openmp_cmd + ' ' + eddy_all_options)
  with open('eddy_output.txt', 'wb') as eddy_output_file:
    eddy_output_file.write((eddy_output.stdout + '\n' + eddy_output.stderr + '\n').encode('utf-8', errors='replace'))
  if app.VERBOSITY > 1:
    app.console('Output of eddy command:')
    sys.stderr.write(eddy_output.stdout + '\n' + eddy_output.stderr + '\n')
  app.cleanup('eddy_in.nii')

  eddy_output_image_path = fsl.find_image('dwi_post_eddy')


  # Check to see whether or not eddy has provided a rotated bvecs file;
  #   if it has, import this into the output image
  bvecs_path = 'dwi_post_eddy.eddy_rotated_bvecs'
  if not os.path.isfile(bvecs_path):
    app.warn('eddy has not provided rotated bvecs file; using original gradient table. Recommend updating FSL eddy to version 5.0.9 or later.')
    bvecs_path = 'bvecs'


  # Run eddy qc tool QUAD if installed and one of -eddyqc_text or -eddyqc_all is specified
  eddyqc_prefix = 'dwi_post_eddy'
  if eddyqc_path:
    if shutil.which('eddy_quad'):

      eddyqc_mask = 'eddy_mask.nii'
      eddyqc_fieldmap = fsl.find_image('field_map') if execute_topup else None
      eddyqc_slspec = 'slspec.txt' if eddy_mporder else None

      # If there was any relevant padding applied, then we want to provide
      #   the comprehensive set of files to EddyQC with that padding removed
      if dwi_post_eddy_crop_option:
        progress = app.ProgressBar('Removing image padding prior to running EddyQC', len(eddyqc_files) + 3)

        for eddy_filename in eddyqc_files:
          if os.path.isfile('dwi_post_eddy.' + eddy_filename):
            if slice_padded and eddy_filename in [ 'eddy_outlier_map', 'eddy_outlier_n_sqr_stdev_map', 'eddy_outlier_n_stdev_map' ]:
              with open('dwi_post_eddy.' + eddy_filename, 'r', encoding='utf-8') as f_eddyfile:
                eddy_data = f_eddyfile.readlines()
              eddy_data_header = eddy_data[0]
              eddy_data = eddy_data[1:]
              for line in eddy_data:
                line = ' '.join(line.strip().split(' ')[:-1])
              with open('dwi_post_eddy_unpad.' + eddy_filename, 'w', encoding='utf-8') as f_eddyfile:
                f_eddyfile.write(eddy_data_header + '\n')
                f_eddyfile.write('\n'.join(eddy_data) + '\n')
            elif eddy_filename.endswith('.nii.gz'):
              run.command('mrconvert dwi_post_eddy.' + eddy_filename + ' dwi_post_eddy_unpad.' + eddy_filename + dwi_post_eddy_crop_option)
            else:
              run.function(os.symlink, 'dwi_post_eddy.' + eddy_filename, 'dwi_post_eddy_unpad.' + eddy_filename)
            app.cleanup('dwi_post_eddy.' + eddy_filename)
          progress.increment()

        if eddy_mporder and slice_padded:
          app.debug('Current slice groups: ' + str(slice_groups))
          app.debug('Slice encoding direction: ' + str(slice_encoding_direction))
          # Remove padded slice from slice_groups, write new slspec
          if sum(slice_encoding_direction) < 0:
            slice_groups = [ [ index-1 for index in group if index ] for group in slice_groups ]
          else:
            slice_groups = [ [ index for index in group if index != dwi_num_slices-1 ] for group in slice_groups ]
          eddyqc_slspec = 'slspec_unpad.txt'
          app.debug('Slice groups after removal: ' + str(slice_groups))
          try:
            # After this removal, slspec should now be a square matrix
            assert all(len(group) == len(slice_groups[0]) for group in slice_groups[1:])
            matrix.save_matrix(eddyqc_slspec, slice_groups, add_to_command_history=False, fmt='%d')
          except AssertionError:
            matrix.save_numeric(eddyqc_slspec, slice_groups, add_to_command_history=False, fmt='%d')
            raise

        run.command('mrconvert eddy_mask.nii eddy_mask_unpad.nii' + dwi_post_eddy_crop_option)
        eddyqc_mask = 'eddy_mask_unpad.nii'
        progress.increment()
        run.command('mrconvert ' + fsl.find_image('field_map') + ' field_map_unpad.nii' + dwi_post_eddy_crop_option)
        eddyqc_fieldmap = 'field_map_unpad.nii'
        progress.increment()
        run.command('mrconvert ' + eddy_output_image_path + ' dwi_post_eddy_unpad.nii.gz' + dwi_post_eddy_crop_option)
        eddyqc_prefix = 'dwi_post_eddy_unpad'
        progress.done()

      eddyqc_options = ' -idx eddy_indices.txt -par eddy_config.txt -b bvals -m ' + eddyqc_mask
      if os.path.isfile(eddyqc_prefix + '.eddy_residuals.nii.gz'):
        eddyqc_options += ' -g ' + bvecs_path
      if execute_topup:
        eddyqc_options += ' -f ' + eddyqc_fieldmap
      if eddy_mporder:
        eddyqc_options += ' -s ' + eddyqc_slspec
      if app.VERBOSITY > 2:
        eddyqc_options += ' -v'
      try:
        run.command('eddy_quad ' + eddyqc_prefix + eddyqc_options)
      except run.MRtrixCmdError as exception:
        with open('eddy_quad_failure_output.txt', 'wb') as eddy_quad_output_file:
          eddy_quad_output_file.write(str(exception).encode('utf-8', errors='replace'))
        app.debug(str(exception))
        app.warn('Error running automated EddyQC tool \'eddy_quad\'; QC data written to "' + eddyqc_path + '" will be files from "eddy" only')
        # Delete the directory if the script only made it partway through
        try:
          shutil.rmtree(eddyqc_prefix + '.qc')
        except OSError:
          pass
    else:
      app.console('Command \'eddy_quad\' not found in PATH; skipping')


  # Have to retain these images until after eddyQC is run
  # If using -eddyqc_all, also write the mask provided to eddy to the output directory;
  #   therefore don't delete it yet here
  if not app.ARGS.eddyqc_all:
    app.cleanup('eddy_mask.nii')
  if execute_topup:
    app.cleanup(fsl.find_image('field_fieldcoef'))


  # Get the axis strides from the input series, so the output image can be modified to match
  stride_option = ' -strides ' + ','.join([str(i) for i in dwi_header.strides()])


  # Determine whether or not volume recombination should be performed
  # This could be either due to use of -rpe_all option, or just due to the data provided with -rpe_header
  # Rather than trying to re-use the code that was used in the case of -rpe_all, run fresh code
  # The phase-encoding scheme needs to be checked also
  volume_matchings = [ dwi_num_volumes ] * dwi_num_volumes
  volume_pairs = [ ]
  app.debug('Commencing gradient direction matching; ' + str(dwi_num_volumes) + ' volumes')
  for index1 in range(dwi_num_volumes):
    if volume_matchings[index1] == dwi_num_volumes: # As yet unpaired
      for index2 in range(index1+1, dwi_num_volumes):
        if volume_matchings[index2] == dwi_num_volumes: # Also as yet unpaired
          # Here, need to check both gradient matching and reversed phase-encode direction
          if not any(dwi_pe_scheme[index1][i] + dwi_pe_scheme[index2][i] for i in range(0,3)) and grads_match(index1, index2):
            volume_matchings[index1] = index2
            volume_matchings[index2] = index1
            volume_pairs.append([index1, index2])
            app.debug('Matched volume ' + str(index1) + ' with ' + str(index2) + '\n' +
                      'Phase encoding: ' + str(dwi_pe_scheme[index1]) + ' ' + str(dwi_pe_scheme[index2]) + '\n' +
                      'Gradients: ' + str(grad[index1]) + ' ' + str(grad[index2]))
            break


  if len(volume_pairs) != int(dwi_num_volumes/2):

    if execute_topup:
      app.cleanup('topup_in.nii')
      app.cleanup(fsl.find_image('field_map'))

    # Convert the resulting volume to the output image, and re-insert the diffusion encoding
    run.command('mrconvert ' + eddy_output_image_path + ' result.mif' + dwi_permvols_posteddy_option + dwi_post_eddy_crop_option + stride_option + ' -fslgrad ' + bvecs_path + ' bvals')
    app.cleanup(eddy_output_image_path)

  else:
    app.console('Detected matching DWI volumes with opposing phase encoding; performing explicit volume recombination')

    # Perform a manual combination of the volumes output by eddy, since LSR is disabled

    # Generate appropriate bvecs / bvals files
    # Particularly if eddy has provided rotated bvecs, since we're combining two volumes into one that
    #   potentially have subject rotation between them (and therefore the sensitisation direction is
    #   not precisely equivalent), the best we can do is take the mean of the two vectors.
    # Manual recombination of volumes needs to take into account the explicit volume matching

    bvecs = matrix.load_matrix(bvecs_path)
    bvecs_combined_transpose = [ ]
    bvals_combined = [ ]

    for pair in volume_pairs:
      bvec_mean = [ 0.5*(bvecs[0][pair[0]] + bvecs[0][pair[1]]),
                    0.5*(bvecs[1][pair[0]] + bvecs[1][pair[1]]),
                    0.5*(bvecs[2][pair[0]] + bvecs[2][pair[1]]) ]
      norm2 = matrix.dot(bvec_mean, bvec_mean)

      # If one diffusion sensitisation gradient direction is reversed with respect to
      #   the other, still want to enable their recombination; but need to explicitly
      #   account for this when averaging the two directions
      if norm2 < 0.5:
        bvec_mean = [ 0.5*(bvecs[0][pair[0]] - bvecs[0][pair[1]]),
                      0.5*(bvecs[1][pair[0]] - bvecs[1][pair[1]]),
                      0.5*(bvecs[2][pair[0]] - bvecs[2][pair[1]]) ]
        norm2 = matrix.dot(bvec_mean, bvec_mean)

      # Occasionally a b=0 volume can have a zero vector
      if norm2:
        factor = 1.0 / math.sqrt(norm2)
        new_vec = [ bvec_mean[0]*factor, bvec_mean[1]*factor, bvec_mean[2]*factor ]
      else:
        new_vec = [ 0.0, 0.0, 0.0 ]
      bvecs_combined_transpose.append(new_vec)
      bvals_combined.append(0.5 * (grad[pair[0]][3] + grad[pair[1]][3]))

    bvecs_combined = matrix.transpose(bvecs_combined_transpose)
    matrix.save_matrix('bvecs_combined', bvecs_combined, add_to_command_history=False)
    matrix.save_vector('bvals_combined', bvals_combined, add_to_command_history=False)

    # Prior to 5.0.8, a bug resulted in the output field map image from topup having an identity transform,
    #   regardless of the transform of the input image
    # Detect this, and manually replace the transform if necessary
    #   (even if this doesn't cause an issue with the subsequent mrcalc command, it may in the future, it's better for
    #   visualising the script intermediate files, and it gives the user a warning about an out-of-date FSL)
    field_map_image = fsl.find_image('field_map')
    field_map_header = image.Header(field_map_image)
    if not image.match('topup_in.nii', field_map_header, up_to_dim=3):
      app.warn('topup output field image has erroneous header; recommend updating FSL to version 5.0.8 or later')
      new_field_map_image = 'field_map_fix.mif'
      run.command('mrtransform ' + field_map_image + ' -replace topup_in.nii ' + new_field_map_image)
      app.cleanup(field_map_image)
      field_map_image = new_field_map_image
    # In FSL 6.0.0, field map image is erroneously constructed with the same number of volumes as the input image,
    #   with all but the first volume containing intensity-scaled duplicates of the uncorrected input images
    # The first volume is however the expected field offset image
    elif len(field_map_header.size()) == 4:
      app.console('Correcting erroneous FSL 6.0.0 field map image output')
      new_field_map_image = 'field_map_fix.mif'
      run.command('mrconvert ' + field_map_image + ' -coord 3 0 -axes 0,1,2 ' + new_field_map_image)
      app.cleanup(field_map_image)
      field_map_image = new_field_map_image
    app.cleanup('topup_in.nii')


    # Derive the weight images
    # Scaling term for field map is identical to the bandwidth provided in the topup config file
    #   (converts Hz to pixel count; that way a simple image gradient can be used to get the Jacobians)
    # Let mrfilter apply the default 1 voxel size gaussian smoothing filter before calculating the field gradient
    #
    #   The jacobian image may be different for any particular volume pair
    #   The appropriate PE directions and total readout times can be acquired from the eddy-style config/index files
    #   eddy_config.txt and eddy_indices.txt
    eddy_config = matrix.load_matrix('eddy_config.txt')
    eddy_indices = matrix.load_vector('eddy_indices.txt', dtype=int)
    app.debug('EDDY config: ' + str(eddy_config))
    app.debug('EDDY indices: ' + str(eddy_indices))

    # This section derives, for each phase encoding configuration present, the 'weight' to be applied
    #   to the image during volume recombination, which is based on the Jacobian of the field in the
    #   phase encoding direction
    for index, config in enumerate(eddy_config):
      pe_axis = [ i for i, e in enumerate(config[0:3]) if e != 0][0]
      sign_multiplier = ' -1.0 -mult' if config[pe_axis] < 0 else ''
      field_derivative_path = 'field_deriv_pe_' + str(index+1) + '.mif'
      run.command('mrcalc ' + field_map_image + ' ' + str(config[3]) + ' -mult' + sign_multiplier + ' - | mrfilter - gradient - | mrconvert - ' + field_derivative_path + ' -coord 3 ' + str(pe_axis) + ' -axes 0,1,2')
      jacobian_path = 'jacobian_' + str(index+1) + '.mif'
      run.command('mrcalc 1.0 ' + field_derivative_path + ' -add 0.0 -max ' + jacobian_path)
      app.cleanup(field_derivative_path)
      run.command('mrcalc ' + jacobian_path + ' ' + jacobian_path + ' -mult weight' + str(index+1) + '.mif')
      app.cleanup(jacobian_path)
    app.cleanup(field_map_image)

    # If eddy provides its main image output in a compressed format, the code block below will need to
    #   uncompress that image independently for every volume pair. Instead, if this is the case, let's
    #   convert it to an uncompressed format before we do anything with it.
    if eddy_output_image_path.endswith('.gz'):
      new_eddy_output_image_path = 'dwi_post_eddy_uncompressed.mif'
      run.command('mrconvert ' + eddy_output_image_path + ' ' + new_eddy_output_image_path)
      app.cleanup(eddy_output_image_path)
      eddy_output_image_path = new_eddy_output_image_path

    # If the DWI volumes were permuted prior to running eddy, then the simplest approach is to permute them
    #   back to their original positions; otherwise, the stored gradient vector directions / phase encode
    #   directions / matched volume pairs are no longer appropriate
    if dwi_permvols_posteddy_option:
      new_eddy_output_image_path = os.path.splitext(eddy_output_image_path)[0] + '_volpermuteundo.mif'
      run.command('mrconvert ' + eddy_output_image_path + dwi_permvols_posteddy_option + ' ' + new_eddy_output_image_path)
      app.cleanup(eddy_output_image_path)
      eddy_output_image_path = new_eddy_output_image_path

    # This section extracts the two volumes corresponding to each reversed phase-encoded volume pair, and
    #   derives a single image volume based on the recombination equation
    combined_image_list = [ ]
    progress = app.ProgressBar('Performing explicit volume recombination', len(volume_pairs))
    for index, volumes in enumerate(volume_pairs):
      pe_indices = [ eddy_indices[i] for i in volumes ]
      run.command('mrconvert ' + eddy_output_image_path + ' volume0.mif -coord 3 ' + str(volumes[0]))
      run.command('mrconvert ' + eddy_output_image_path + ' volume1.mif -coord 3 ' + str(volumes[1]))
      # Volume recombination equation described in Skare and Bammer 2010
      combined_image_path = 'combined' + str(index) + '.mif'
      run.command('mrcalc volume0.mif weight' + str(pe_indices[0]) + '.mif -mult volume1.mif weight' + str(pe_indices[1]) + '.mif -mult -add weight' + str(pe_indices[0]) + '.mif weight' + str(pe_indices[1]) + '.mif -add -divide 0.0 -max ' + combined_image_path)
      combined_image_list.append(combined_image_path)
      run.function(os.remove, 'volume0.mif')
      run.function(os.remove, 'volume1.mif')
      progress.increment()
    progress.done()

    app.cleanup(eddy_output_image_path)
    for index in range(0, len(eddy_config)):
      app.cleanup('weight' + str(index+1) + '.mif')

    # Finally the recombined volumes must be concatenated to produce the resulting image series
    combine_command = ['mrcat', combined_image_list, '-', '-axis', '3', '|', \
                       'mrconvert', '-', 'result.mif', '-fslgrad', 'bvecs_combined', 'bvals_combined']
    if dwi_post_eddy_crop_option:
      combine_command.extend(dwi_post_eddy_crop_option.strip().split(' '))
    combine_command.extend(stride_option.strip().split(' '))
    run.command(combine_command)
    app.cleanup(combined_image_list)


  # Grab any relevant files that eddy has created, and copy them to the requested directory
  if eddyqc_path:
    if app.FORCE_OVERWRITE and os.path.exists(eddyqc_path) and not os.path.isdir(eddyqc_path):
      run.function(os.remove, eddyqc_path)
    if not os.path.exists(eddyqc_path):
      run.function(os.makedirs, eddyqc_path)
    for filename in eddyqc_files:
      if os.path.exists(eddyqc_prefix + '.' + filename):
        # If this is an image, and axis padding was applied, want to undo the padding
        if filename.endswith('.nii.gz') and dwi_post_eddy_crop_option:
          run.command('mrconvert ' + eddyqc_prefix + '.' + filename + ' ' + shlex.quote(os.path.join(eddyqc_path, filename)) + dwi_post_eddy_crop_option, force=app.FORCE_OVERWRITE)
        else:
          run.function(shutil.copy, eddyqc_prefix + '.' + filename, os.path.join(eddyqc_path, filename))
    # Also grab any files generated by the eddy qc tool QUAD
    if os.path.isdir(eddyqc_prefix + '.qc'):
      if app.FORCE_OVERWRITE and os.path.exists(os.path.join(eddyqc_path, 'quad')):
        run.function(shutil.rmtree, os.path.join(eddyqc_path, 'quad'))
      run.function(shutil.copytree, eddyqc_prefix + '.qc', os.path.join(eddyqc_path, 'quad'))
    # Also grab the brain mask that was provided to eddy if -eddyqc_all was specified
    if app.ARGS.eddyqc_all:
      if dwi_post_eddy_crop_option:
        run.command('mrconvert eddy_mask.nii ' + shlex.quote(os.path.join(eddyqc_path, 'eddy_mask.nii')) + dwi_post_eddy_crop_option, force=app.FORCE_OVERWRITE)
      else:
        run.function(shutil.copy, 'eddy_mask.nii', os.path.join(eddyqc_path, 'eddy_mask.nii'))
      app.cleanup('eddy_mask.nii')




  keys_to_remove = [ 'MultibandAccelerationFactor', 'SliceEncodingDirection', 'SliceTiming' ]
  # These keys are still relevant for the output data if no EPI distortion correction was performed
  if execute_applytopup:
    keys_to_remove.extend([ 'PhaseEncodingDirection', 'TotalReadoutTime', 'pe_scheme' ])
  # Get the header key-value entries from the input DWI, remove those we don't wish to keep, and
  #   export the result to a new JSON file so that they can be inserted into the output header
  with open('dwi.json', 'r', encoding='utf-8') as input_json_file:
    keyval = json.load(input_json_file)
  for key in keys_to_remove:
    keyval.pop(key, None)
  # Make sure to use the revised diffusion gradient table rather than that of the input;
  #  incorporates motion correction, and possibly also the explicit volume recombination
  keyval['dw_scheme'] = image.Header('result.mif').keyval()['dw_scheme']
  # 'Stash' the phase encoding scheme of the original uncorrected DWIs, since it still
  #   may be useful information at some point in the future but is no longer relevant
  #   for e.g. tracking for different volumes, or performing any geometric corrections
  if execute_applytopup:
    keyval['prior_pe_scheme'] = dwi_manual_pe_scheme if dwi_manual_pe_scheme else dwi_pe_scheme
  with open('output.json', 'w', encoding='utf-8') as output_json_file:
    json.dump(keyval, output_json_file)


  # Finish!
  run.command('mrconvert result.mif ' + path.from_user(app.ARGS.output) + grad_export_option, mrconvert_keyval='output.json', force=app.FORCE_OVERWRITE)

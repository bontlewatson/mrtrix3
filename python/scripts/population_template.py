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

# Generates an unbiased group-average template via image registration of images to a midway space.

import json, math, os, re, shlex, shutil, sys

DEFAULT_RIGID_SCALES  = [0.3,0.4,0.6,0.8,1.0,1.0]
DEFAULT_RIGID_LMAX    = [2,2,2,4,4,4]
DEFAULT_AFFINE_SCALES = [0.3,0.4,0.6,0.8,1.0,1.0]
DEFAULT_AFFINE_LMAX   = [2,2,2,4,4,4]

DEFAULT_NL_SCALES = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
DEFAULT_NL_NITER  = [  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5]
DEFAULT_NL_LMAX   = [  2,  2,  2,  2,  2,  2,  2,  2,  4,  4,  4,  4,  4,  4,  4,  4]

REGISTRATION_MODES = ['rigid', 'affine', 'nonlinear', 'rigid_affine', 'rigid_nonlinear', 'affine_nonlinear', 'rigid_affine_nonlinear']

AGGREGATION_MODES = ['mean', 'median']

IMAGEEXT = ['mif', 'nii', 'mih', 'mgh', 'mgz', 'img', 'hdr']

def usage(cmdline): #pylint: disable=unused-variable
  cmdline.set_author('David Raffelt (david.raffelt@florey.edu.au) & Max Pietsch (maximilian.pietsch@kcl.ac.uk) & Thijs Dhollander (thijs.dhollander@gmail.com)')

  cmdline.set_synopsis('Generates an unbiased group-average template from a series of images')
  cmdline.add_description('First a template is optimised with linear registration (rigid and/or affine, both by default), then non-linear registration is used to optimise the template further.')
  cmdline.add_argument('input_dir', nargs='+', help='Directory containing all input images of a given contrast')
  cmdline.add_argument('template', help='Output template image')

  cmdline.add_example_usage('Multi-contrast registration',
                            'population_template input_WM_ODFs/ output_WM_template.mif input_GM_ODFs/ output_GM_template.mif',
                            'When performing multi-contrast registration, the input directory and corresponding output template '
                            'image for a given contrast are to be provided as a pair, '
                            'with the pairs corresponding to different contrasts provided sequentially.')

  options = cmdline.add_argument_group('Multi-contrast options')
  options.add_argument('-mc_weight_initial_alignment', help='Weight contribution of each contrast to the initial alignment. Comma separated, default: 1.0')
  options.add_argument('-mc_weight_rigid', help='Weight contribution of each contrast to the objective of rigid registration. Comma separated, default: 1.0')
  options.add_argument('-mc_weight_affine', help='Weight contribution of each contrast to the objective of affine registration. Comma separated, default: 1.0')
  options.add_argument('-mc_weight_nl', help='Weight contribution of each contrast to the objective of nonlinear registration. Comma separated, default: 1.0')

  linoptions = cmdline.add_argument_group('Options for the linear registration')
  linoptions.add_argument('-linear_no_pause', action='store_true', help='Do not pause the script if a linear registration seems implausible')
  linoptions.add_argument('-linear_no_drift_correction', action='store_true', help='Deactivate correction of template appearance (scale and shear) over iterations')
  linoptions.add_argument('-linear_estimator', help='Specify estimator for intensity difference metric. Valid choices are: l1 (least absolute: |x|), l2 (ordinary least squares), lp (least powers: |x|^1.2), Default: None (no robust estimator used)')
  linoptions.add_argument('-rigid_scale', help='Specify the multi-resolution pyramid used to build the rigid template, in the form of a list of scale factors (default: %s). This and affine_scale implicitly  define the number of template levels' % ','.join([str(x) for x in DEFAULT_RIGID_SCALES]))
  linoptions.add_argument('-rigid_lmax', help='Specify the lmax used for rigid registration for each scale factor, in the form of a list of integers (default: %s). The list must be the same length as the linear_scale factor list' % ','.join([str(x) for x in DEFAULT_RIGID_LMAX]))
  linoptions.add_argument('-rigid_niter', help='Specify the number of registration iterations used within each level before updating the template, in the form of a list of integers (default:50 for each scale). This must be a single number or a list of same length as the linear_scale factor list')
  linoptions.add_argument('-affine_scale', help='Specify the multi-resolution pyramid used to build the affine template, in the form of a list of scale factors (default: %s). This and rigid_scale implicitly define the number of template levels' % ','.join([str(x) for x in DEFAULT_AFFINE_SCALES]))
  linoptions.add_argument('-affine_lmax', help='Specify the lmax used for affine registration for each scale factor, in the form of a list of integers (default: %s). The list must be the same length as the linear_scale factor list' % ','.join([str(x) for x in DEFAULT_AFFINE_LMAX]))
  linoptions.add_argument('-affine_niter', help='Specify the number of registration iterations used within each level before updating the template, in the form of a list of integers (default:500 for each scale). This must be a single number or a list of same length as the linear_scale factor list')

  nloptions = cmdline.add_argument_group('Options for the non-linear registration')
  nloptions.add_argument('-nl_scale', help='Specify the multi-resolution pyramid used to build the non-linear template, in the form of a list of scale factors (default: %s). This implicitly defines the number of template levels' % ','.join([str(x) for x in DEFAULT_NL_SCALES]))
  nloptions.add_argument('-nl_lmax', help='Specify the lmax used for non-linear registration for each scale factor, in the form of a list of integers (default: %s). The list must be the same length as the nl_scale factor list' % ','.join([str(x) for x in DEFAULT_NL_LMAX]))
  nloptions.add_argument('-nl_niter', help='Specify the number of registration iterations used within each level before updating the template, in the form of a list of integers (default: %s). The list must be the same length as the nl_scale factor list' % ','.join([str(x) for x in DEFAULT_NL_NITER]))
  nloptions.add_argument('-nl_update_smooth', default='2.0', help='Regularise the gradient update field with Gaussian smoothing (standard deviation in voxel units, Default 2.0 x voxel_size)')
  nloptions.add_argument('-nl_disp_smooth', default='1.0', help='Regularise the displacement field with Gaussian smoothing (standard deviation in voxel units, Default 1.0 x voxel_size)')
  nloptions.add_argument('-nl_grad_step', default='0.5', help='The gradient step size for non-linear registration (Default: 0.5)')

  options = cmdline.add_argument_group('Input, output and general options')
  options.add_argument('-type', help='Specify the types of registration stages to perform. Options are "rigid" (perform rigid registration only which might be useful for intra-subject registration in longitudinal analysis), "affine" (perform affine registration) and "nonlinear" as well as cominations of registration types: %s. Default: rigid_affine_nonlinear' % ', '.join('"' + x + '"' for x in REGISTRATION_MODES if "_" in x), default='rigid_affine_nonlinear')
  options.add_argument('-voxel_size', help='Define the template voxel size in mm. Use either a single value for isotropic voxels or 3 comma separated values.')
  options.add_argument('-initial_alignment', default='mass', help='Method of alignment to form the initial template. Options are "mass" (default), "robust_mass" (requires masks), "geometric" and "none".')
  options.add_argument('-mask_dir', help='Optionally input a set of masks inside a single directory, one per input image (with the same file name prefix). Using masks will speed up registration significantly. Note that masks are used for registration, not for aggregation. To exclude areas from aggregation, NaN-mask your input images.')
  options.add_argument('-warp_dir', help='Output a directory containing warps from each input to the template. If the folder does not exist it will be created')
  options.add_argument('-transformed_dir', help='Output a directory containing the input images transformed to the template. If the folder does not exist it will be created. For multi-contrast registration, provide comma separated list of directories.')
  options.add_argument('-linear_transformations_dir', help='Output a directory containing the linear transformations used to generate the template. If the folder does not exist it will be created')
  options.add_argument('-template_mask', help='Output a template mask. Only works if -mask_dir has been input. The template mask is computed as the intersection of all subject masks in template space.')
  options.add_argument('-noreorientation', action='store_true', help='Turn off FOD reorientation in mrregister. Reorientation is on by default if the number of volumes in the 4th dimension corresponds to the number of coefficients in an antipodally symmetric spherical harmonic series (i.e. 6, 15, 28, 45, 66 etc)')
  options.add_argument('-leave_one_out', help='Register each input image to a template that does not contain that image. Valid choices: 0, 1, auto. (Default: auto (true if n_subjects larger than 2 and smaller than 15)) ')
  options.add_argument('-aggregate', help='Measure used to aggregate information from transformed images to the template image. Valid choices: %s. Default: mean' % ', '.join(AGGREGATION_MODES))
  options.add_argument('-aggregation_weights', help='Comma separated file containing weights used for weighted image aggregation. Each row must contain the identifiers of the input image and its weight. Note that this weighs intensity values not transformations (shape).')
  options.add_argument('-nanmask', action='store_true', help='Optionally apply masks to (transformed) input images using NaN values to specify include areas for registration and aggregation. Only works if -mask_dir has been input.')
  options.add_argument('-copy_input', action='store_true', help='Copy input images and masks into local scratch directory.')
  options.add_argument('-delete_temporary_files', action='store_true', help='Delete temporary files from scratch directory during template creation.')

# ENH: add option to initialise warps / transformations



def abspath(arg, *args):
  return os.path.abspath(os.path.join(arg, *args))


def relpath(arg, *args):
  from mrtrix3 import app #pylint: disable=no-name-in-module, import-outside-toplevel
  return os.path.relpath(os.path.join(arg, *args), app.WORKING_DIR)


def copy(src, dst, follow_symlinks=True):
  """Copy data but do not set mode bits. Return the file's destination.

  mimics shutil.copy but without setting mode bits as shutil.copymode can fail on exotic mounts
  (observed on cifs with file_mode=0777).
  """
  if os.path.isdir(dst):
    dst = os.path.join(dst, os.path.basename(src))
  if sys.version_info[0] > 2:
    shutil.copyfile(src, dst, follow_symlinks=follow_symlinks)   # pylint: disable=unexpected-keyword-arg
  else:
    shutil.copyfile(src, dst)
  return dst


def check_linear_transformation(transformation, cmd, max_scaling=0.5, max_shear=0.2, max_rot=None, pause_on_warn=True):
  from mrtrix3 import app, run, utils #pylint: disable=no-name-in-module, import-outside-toplevel
  if max_rot is None:
    max_rot = 2 * math.pi

  good = True
  run.command('transformcalc ' + transformation + ' decompose ' + transformation + 'decomp')
  if not os.path.isfile(transformation + 'decomp'):  # does not exist if run with -continue option
    app.console(transformation + 'decomp not found. skipping check')
    return True
  data = utils.load_keyval(transformation + 'decomp')
  run.function(os.remove, transformation + 'decomp')
  scaling = [float(value) for value in data['scaling']]
  if any(a < 0 for a in scaling) or any(a > (1 + max_scaling) for a in scaling) or any(
      a < (1 - max_scaling) for a in scaling):
    app.warn("large scaling (" + str(scaling) + ") in " + transformation)
    good = False
  shear = [float(value) for value in data['shear']]
  if any(abs(a) > max_shear for a in shear):
    app.warn("large shear (" + str(shear) + ") in " + transformation)
    good = False
  rot_angle = float(data['angle_axis'][0])
  if abs(rot_angle) > max_rot:
    app.warn("large rotation (" + str(rot_angle) + ") in " + transformation)
    good = False

  if not good:
    newcmd = []
    what = ''
    init_rotation_found = False
    skip = 0
    for element in cmd.split():
      if skip:
        skip -= 1
        continue
      if '_init_rotation' in element:
        init_rotation_found = True
      if '_init_matrix' in element:
        skip = 1
        continue
      if 'affine_scale' in element:
        assert what != 'rigid'
        what = 'affine'
      elif 'rigid_scale' in element:
        assert what != 'affine'
        what = 'rigid'
      newcmd.append(element)
    newcmd = " ".join(newcmd)
    if not init_rotation_found:
      app.console("replacing the transformation obtained with:")
      app.console(cmd)
      if what:
        newcmd += ' -' + what + '_init_translation mass -' + what + '_init_rotation search'
      app.console("by the one obtained with:")
      app.console(newcmd)
      run.command(newcmd, force=True)
      return check_linear_transformation(transformation, newcmd, max_scaling, max_shear, max_rot, pause_on_warn=pause_on_warn)
    if pause_on_warn:
      app.warn("you might want to manually repeat mrregister with different parameters and overwrite the transformation file: \n%s" % transformation)
      app.console('The command that failed the test was: \n' + cmd)
      app.console('Working directory: \n' + os.getcwd())
      input("press enter to continue population_template")
  return good


def aggregate(inputs, output, contrast_idx, mode, force=True):
  from mrtrix3 import MRtrixError, run  # pylint: disable=no-name-in-module, import-outside-toplevel

  images = [inp.ims_transformed[contrast_idx] for inp in inputs]
  if mode == 'mean':
    run.command(['mrmath', images, 'mean', '-keep_unary_axes', output], force=force)
  elif mode == 'median':
    run.command(['mrmath', images, 'median', '-keep_unary_axes', output], force=force)
  elif mode == 'weighted_mean':
    weights = [inp.aggregation_weight for inp in inputs]
    assert not any(w is None for w in weights), weights
    wsum = sum(float(w) for w in weights)
    cmd = ['mrcalc']
    if wsum <= 0:
      raise MRtrixError("the sum of aggregetion weights has to be positive")
    for weight, image in zip(weights, images):
      if float(weight) != 0:
        cmd += [image, weight, '-mult'] + (['-add'] if len(cmd) > 1 else [])
    cmd += ['%.16f' % wsum, '-div', output]
    run.command(cmd, force=force)
  else:
    raise MRtrixError("aggregation mode %s not understood" % mode)


def inplace_nan_mask(images, masks):
  from mrtrix3 import run  # pylint: disable=no-name-in-module, import-outside-toplevel
  assert len(images) == len(masks), (len(images), len(masks))
  for image, mask in zip(images, masks):
    target_dir = os.path.split(image)[0]
    masked = os.path.join(target_dir, '__' + os.path.split(image)[1])
    run.command("mrcalc " + mask + " " + image + " nan -if " + masked, force=True)
    run.function(shutil.move, masked, image)


def calculate_isfinite(inputs, contrasts):
  from mrtrix3 import run, path  # pylint: disable=no-name-in-module, import-outside-toplevel
  agg_weights = [float(inp.aggregation_weight) for inp in inputs if inp.aggregation_weight is not None]
  for cid in range(contrasts.n_contrasts):
    for inp in inputs:
      if contrasts.n_volumes[cid] > 0:
        cmd = 'mrconvert ' + inp.ims_transformed[cid] + ' -coord 3 0 - | mrcalc - -finite'
      else:
        cmd = 'mrcalc ' + inp.ims_transformed[cid] + ' -finite'
      if inp.aggregation_weight:
        cmd += ' %s -mult ' % inp.aggregation_weight
      cmd += ' isfinite%s/%s.mif' % (contrasts.suff[cid], inp.uid)
      run.command(cmd, force=True)
  for cid in range(contrasts.n_contrasts):
    cmd = ['mrmath', path.all_in_dir('isfinite%s' % contrasts.suff[cid]), 'sum']
    if agg_weights:
      agg_weight_norm = str(float(len(agg_weights)) / sum(agg_weights))
      cmd += ['-', '|', 'mrcalc', '-', agg_weight_norm, '-mult']
    run.command(cmd + [contrasts.isfinite_count[cid]], force=True)


def get_common_postfix(file_list):
  return os.path.commonprefix([i[::-1] for i in file_list])[::-1]


def get_common_prefix(file_list):
  return os.path.commonprefix(file_list)


class Contrasts:
  """
      Class that parses arguments and holds information specific to each image contrast

      Attributes
      ----------
      suff: list of str
        identifiers used for contrast-specific filenames and folders ['_c0', '_c1', ...]

      names: list of str
        derived from constrast-specific input folder

      templates_out: list of str
        full path to output templates

      templates: list of str
        holds current template names during registration

      n_volumes: list of int
        number of volumes in each contrast

      fod_reorientation: list of bool
        whether to perform FOD reorientation with mrtransform

      isfinite_count: list of str
        filenames of images holding (weighted) number of finite-valued voxels across all images

      mc_weight_<mode>: list of str
        contrast-specific weight used during initialisation / registration

      <mode>_weight_option: list of str
        weight option to be passed to mrregister, <mode> = {'initial_alignment', 'rigid', 'affine', 'nl'}

      n_contrasts: int

      """

  def __init__(self):
    from mrtrix3 import MRtrixError, path, app  # pylint: disable=no-name-in-module, import-outside-toplevel

    n_contrasts = len(app.ARGS.input_dir)

    self.suff = ["_c" + c for c in map(str, range(n_contrasts))]
    self.names = [os.path.relpath(f, os.path.commonprefix(app.ARGS.input_dir)) for f in app.ARGS.input_dir]

    self.templates_out = [path.from_user(t, True) for t in app.ARGS.template]

    self.mc_weight_initial_alignment = [None for _ in range(self.n_contrasts)]
    self.mc_weight_rigid = [None for _ in range(self.n_contrasts)]
    self.mc_weight_affine = [None for _ in range(self.n_contrasts)]
    self.mc_weight_nl = [None for _ in range(self.n_contrasts)]
    self.initial_alignment_weight_option = [None for _ in range(self.n_contrasts)]
    self.rigid_weight_option = [None for _ in range(self.n_contrasts)]
    self.affine_weight_option = [None for _ in range(self.n_contrasts)]
    self.nl_weight_option = [None for _ in range(self.n_contrasts)]

    self.isfinite_count = ['isfinite' + c + '.mif' for c in self.suff]
    self.templates = [None for _ in range(self.n_contrasts)]
    self.n_volumes = [None for _ in range(self.n_contrasts)]
    self.fod_reorientation = [None for _ in range(self.n_contrasts)]


    for mode in ['initial_alignment', 'rigid', 'affine', 'nl']:
      opt = app.ARGS.__dict__.get('mc_weight_' + mode, None)
      if opt:
        if n_contrasts == 1:
          raise MRtrixError('mc_weight_' + mode+' requires multiple input contrasts')
        opt = opt.split(',')
        if len(opt) != n_contrasts:
          raise MRtrixError('mc_weight_' + mode+' needs to be defined for each contrast')
      else:
        opt = ["1"] * n_contrasts
      self.__dict__['mc_weight_%s' % mode] = opt
      self.__dict__['%s_weight_option' % mode] = ' -mc_weights '+','.join(opt)+' ' if n_contrasts > 1 else ''

    if len(self.templates_out) != n_contrasts:
      raise MRtrixError('number of templates (%i) does not match number of input directories (%i)' %
                        (len(self.templates_out), n_contrasts))

  @property
  def n_contrasts(self):
    return len(self.suff)

  def __repr__(self, *args, **kwargs):
    text = ''
    for cid in range(self.n_contrasts):
      text += '\tcontrast: %s, template: %s, suffix: %s\n' % (self.names[cid], self.templates_out[cid], self.suff[cid])
    return text


class Input:
  """
      Class that holds input information specific to a single image (multiple contrasts)

      Attributes
      ----------
      uid: str
        unique identifier for these input image(s), does not contain spaces

      ims_path: list of str
        full path to input images, shell quoted OR paths to cached file if cache_local was called

      msk_path: str
        full path to input mask, shell quoted OR path to cached file if cache_local was called

      ims_filenames : list of str
        for each contrast the input file paths stripped of their respective directories. Used for final output only.

      msk_filename: str
        as ims_filenames

      ims_transformed: list of str
        input_transformed<contrast identifier>/<uid>.mif

      msk_transformed: list of str
        mask_transformed/<uid>.mif

      aggregation_weight: float
        weights used in image aggregation that forms the template. Has to be normalised across inputs.

      _im_directories : list of str
        full path to user-provided input directories containing the input images, one for each contrast

      _msk_directory: str
        full path to user-provided mask directory

      _local_ims: list of str
        path to cached input images

      _local_msk: str
        path to cached input mask

      Methods
      -------
      cache_local()
        copy files into folders in current working directory. modifies _local_ims and  _local_msk

      """
  def __init__(self, uid, filenames, directories, contrasts, mask_filename='', mask_directory=''):
    self.contrasts = contrasts

    self.uid = uid
    assert self.uid, "UID empty"
    assert self.uid.count(' ') == 0, 'UID "%s" contains whitespace' % self.uid

    assert len(directories) == len(filenames)
    self.ims_filenames = filenames
    self._im_directories = directories

    self.msk_filename = mask_filename
    self._msk_directory = mask_directory

    n_contrasts = len(contrasts)

    self.ims_transformed = [os.path.join('input_transformed'+contrasts[cid], uid + '.mif') for cid in range(n_contrasts)]
    self.msk_transformed = os.path.join('mask_transformed', uid + '.mif')

    self.aggregation_weight = None

    self._local_ims = []
    self._local_msk = None

  def __repr__(self, *args, **kwargs):
    text = '\nInput ['
    for key in sorted([k for k in self.__dict__ if not k.startswith('_')]):
      text += '\n\t' + str(key) + ': ' + str(self.__dict__[key])
    text += '\n]'
    return text

  def info(self):
    message = ['input: ' + self.uid]
    if self.aggregation_weight:
      message += ['agg weight: ' + self.aggregation_weight]
    for csuff, fname in zip(self.contrasts, self.ims_filenames):
      message += [((csuff + ': ') if csuff else '') + '"' + fname + '"']
    if self.msk_filename:
      message += ['mask: ' + self.msk_filename]
    return ', '.join(message)

  def cache_local(self):
    from mrtrix3 import run, path  # pylint: disable=no-name-in-module, import-outside-toplevel
    contrasts = self.contrasts
    for cid, csuff in enumerate(contrasts):
      if not os.path.isdir('input' + csuff):
        path.make_dir('input' + csuff)
      run.command('mrconvert ' + self.ims_path[cid] + ' ' + os.path.join('input' + csuff, self.uid + '.mif'))
    self._local_ims = [os.path.join('input' + csuff, self.uid + '.mif') for csuff in contrasts]
    if self.msk_filename:
      if not os.path.isdir('mask'):
        path.make_dir('mask')
      run.command('mrconvert ' + self.msk_path + ' ' + os.path.join('mask', self.uid + '.mif'))
      self._local_msk = os.path.join('mask', self.uid + '.mif')

  def get_ims_path(self, quoted=True):
    """ return path to input images """
    from mrtrix3 import path  # pylint: disable=no-name-in-module, import-outside-toplevel
    if self._local_ims:
      return self._local_ims
    return [path.from_user(abspath(d, f), quoted) for d, f in zip(self._im_directories, self.ims_filenames)]
  ims_path = property(get_ims_path)

  def get_msk_path(self, quoted=True):
    """ return path to input mask """
    from mrtrix3 import path  # pylint: disable=no-name-in-module, import-outside-toplevel
    if self._local_msk:
      return self._local_msk
    return path.from_user(os.path.join(self._msk_directory, self.msk_filename), quoted) if self.msk_filename else None
  msk_path = property(get_msk_path)


def parse_input_files(in_files, mask_files, contrasts, f_agg_weight=None, whitespace_repl='_'):
  """
    matches input images across contrasts and pair them with masks.
    extracts unique identifiers from mask and image filenames by stripping common pre and postfix (per contrast and for masks)
    unique identifiers contain ASCII letters, numbers and '_' but no whitespace which is replaced by whitespace_repl

    in_files: list of lists
      the inner list holds filenames specific to a contrast

    mask_files:
      can be empty

    returns list of Input

    checks: 3d_nonunity
    TODO check if no common grid & trafo across contrasts (only relevant for robust init?)

  """
  from mrtrix3 import MRtrixError, app, path, image  # pylint: disable=no-name-in-module, import-outside-toplevel
  contrasts = contrasts.suff
  inputs = []
  def paths_to_file_uids(paths, prefix, postfix):
    """ strip pre and postfix from filename, replace whitespace characters """
    uid_path = {}
    uids = []
    for path in paths:
      uid = re.sub(re.escape(postfix)+'$', '', re.sub('^'+re.escape(prefix), '', os.path.split(path)[1]))
      uid = re.sub(r'\s+', whitespace_repl, uid)
      if not uid:
        raise MRtrixError('No uniquely identifiable part of filename "' + path + '" '
                          'after prefix and postfix substitution '
                          'with prefix "' + prefix + '" and postfix "' + postfix + '"')
      app.debug('UID mapping: "' + path + '" --> "' + uid + '"')
      if uid in uid_path:
        raise MRtrixError('unique file identifier is not unique: "' + uid + '" mapped to "' + path + '" and "' + uid_path[uid] +'"')
      uid_path[uid] = path
      uids.append(uid)
    return uids

  # mask uids
  mask_uids = []
  if mask_files:
    mask_common_postfix = get_common_postfix(mask_files)
    if not mask_common_postfix:
      raise MRtrixError('mask filenames do not have a common postfix')
    mask_common_prefix = get_common_prefix([os.path.split(m)[1] for m in mask_files])
    mask_uids = paths_to_file_uids(mask_files, mask_common_prefix, mask_common_postfix)
    if app.VERBOSITY > 1:
      app.console('mask uids:' + str(mask_uids))

  # images uids
  common_postfix = [get_common_postfix(files) for files in in_files]
  common_prefix = [get_common_prefix(files) for files in in_files]
  # xcontrast_xsubject_pre_postfix: prefix and postfix of the common part across contrasts and subjects,
  # without image extensions and leading or trailing '_' or '-'
  xcontrast_xsubject_pre_postfix = [get_common_postfix(common_prefix).lstrip('_-'),
                                    get_common_prefix([re.sub('.('+'|'.join(IMAGEEXT)+')(.gz)?$', '', pfix).rstrip('_-') for pfix in common_postfix])]
  if app.VERBOSITY > 1:
    app.console("common_postfix: " + str(common_postfix))
    app.console("common_prefix: " + str(common_prefix))
    app.console("xcontrast_xsubject_pre_postfix: " + str(xcontrast_xsubject_pre_postfix))
  for ipostfix, postfix in enumerate(common_postfix):
    if not postfix:
      raise MRtrixError('image filenames do not have a common postfix:\n' + '\n'.join(in_files[ipostfix]))

  c_uids = []
  for cid, files in enumerate(in_files):
    c_uids.append(paths_to_file_uids(files, common_prefix[cid], common_postfix[cid]))

  if app.VERBOSITY > 1:
    app.console('uids by contrast:' + str(c_uids))

  # join images and masks
  for ifile, fname in enumerate(in_files[0]):
    uid = c_uids[0][ifile]
    fnames = [fname]
    dirs = [abspath(path.from_user(app.ARGS.input_dir[0], False))]
    if len(contrasts) > 1:
      for cid in range(1, len(contrasts)):
        dirs.append(abspath(path.from_user(app.ARGS.input_dir[cid], False)))
        image.check_3d_nonunity(os.path.join(dirs[cid], in_files[cid][ifile]))
        if uid != c_uids[cid][ifile]:
          raise MRtrixError('no matching image was found for image %s and contrasts %s and %s.' % (fname, dirs[0], dirs[cid]))
        fnames.append(in_files[cid][ifile])

    if mask_files:
      if uid not in mask_uids:
        raise MRtrixError('no matching mask image was found for input image ' + fname + ' with uid "'+uid+'". '
                          'Mask uid candidates: ' + ', '.join(['"%s"' % m for m in mask_uids]))
      index = mask_uids.index(uid)
      # uid, filenames, directories, contrasts, mask_filename = '', mask_directory = '', agg_weight = None
      inputs.append(Input(uid, fnames, dirs, contrasts,
                          mask_filename=mask_files[index], mask_directory=abspath(path.from_user(app.ARGS.mask_dir, False))))
    else:
      inputs.append(Input(uid, fnames, dirs, contrasts))

  # parse aggregation weights and match to inputs
  if f_agg_weight:
    import csv  # pylint: disable=import-outside-toplevel
    try:
      with open(f_agg_weight, 'r', encoding='utf-8') as fweights:
        agg_weights = dict((row[0].lstrip().rstrip(), row[1]) for row in csv.reader(fweights, delimiter=',', quotechar='#'))
    except UnicodeDecodeError:
      with open(f_agg_weight, 'r', encoding='utf-8') as fweights:
        reader = csv.reader(fweights.read().decode('utf-8', errors='replace'), delimiter=',', quotechar='#')
        agg_weights = dict((row[0].lstrip().rstrip(), row[1]) for row in reader)
    pref = '^' + re.escape(get_common_prefix(list(agg_weights.keys())))
    suff = re.escape(get_common_postfix(list(agg_weights.keys()))) + '$'
    for key in agg_weights.keys():
      agg_weights[re.sub(suff, '', re.sub(pref, '', key))] = agg_weights.pop(key).strip()

    for inp in inputs:
      if inp.uid not in agg_weights:
        raise MRtrixError('aggregation weight not found for %s' % inp.uid)
      inp.aggregation_weight = agg_weights[inp.uid]
    app.console('Using aggregation weights ' + f_agg_weight)
    weights = [float(inp.aggregation_weight) for inp in inputs if inp.aggregation_weight is not None]
    if sum(weights) <= 0:
      raise MRtrixError('Sum of aggregation weights is not positive: ' + str(weights))
    if any(w < 0 for w in weights):
      app.warn('Negative aggregation weights: ' + str(weights))

  return inputs, xcontrast_xsubject_pre_postfix


def execute(): #pylint: disable=unused-variable
  from mrtrix3 import MRtrixError, app, image, matrix, path, run, EXE_LIST #pylint: disable=no-name-in-module, import-outside-toplevel

  expected_commands = ['mrgrid', 'mrregister', 'mrtransform', 'mraverageheader', 'mrconvert', 'mrmath', 'transformcalc', 'mrfilter']
  for cmd in expected_commands:
    if cmd not in EXE_LIST :
      raise MRtrixError("Could not find " + cmd + " in bin/. Binary commands not compiled?")

  if not app.ARGS.type in REGISTRATION_MODES:
    raise MRtrixError("registration type must be one of %s. provided: %s" % (str(REGISTRATION_MODES), app.ARGS.type))
  dorigid     = "rigid"     in app.ARGS.type
  doaffine    = "affine"    in app.ARGS.type
  dolinear    = dorigid or doaffine
  dononlinear = "nonlinear" in app.ARGS.type
  assert (dorigid + doaffine + dononlinear >= 1), "FIXME: registration type not valid"


  input_output = app.ARGS.input_dir + [app.ARGS.template]
  n_contrasts = len(input_output) // 2
  if len(input_output) != 2 * n_contrasts:
    raise MRtrixError('expected two arguments per contrast, received %i: %s' % (len(input_output), ', '.join(input_output)))
  if n_contrasts > 1:
    app.console('Generating population template using multi-contrast registration')

  # reorder arguments for multi-contrast registration as after command line parsing app.ARGS.input_dir holds all but one argument
  app.ARGS.input_dir = []
  app.ARGS.template = []
  for i_contrast in range(n_contrasts):
    inargs = (input_output[i_contrast*2], input_output[i_contrast*2+1])
    if not os.path.isdir(inargs[0]):
      raise MRtrixError('input directory %s not found' % inargs[0])
    app.ARGS.input_dir.append(relpath(inargs[0]))
    app.ARGS.template.append(relpath(inargs[1]))

  cns = Contrasts()
  app.debug(str(cns))

  in_files = [sorted(path.all_in_dir(input_dir, dir_path=False)) for input_dir in app.ARGS.input_dir]
  if len(in_files[0]) <= 1:
    raise MRtrixError('Not enough images found in input directory ' + app.ARGS.input_dir[0] +
                      '. More than one image is needed to generate a population template')
  if n_contrasts > 1:
    for cid in range(1, n_contrasts):
      if len(in_files[cid]) != len(in_files[0]):
        raise MRtrixError('Found %i images in input directory %s ' % (len(app.ARGS.input_dir[0]), app.ARGS.input_dir[0]) +
                          'but %i input images in %s.' % (len(app.ARGS.input_dir[cid]), app.ARGS.input_dir[cid]))
  else:
    app.console('Generating a population-average template from ' + str(len(in_files[0])) + ' input images')
    if n_contrasts > 1:
      app.console('using ' + str(len(in_files)) + ' contrasts for each input image')

  voxel_size = None
  if app.ARGS.voxel_size:
    voxel_size = app.ARGS.voxel_size.split(',')
    if len(voxel_size) == 1:
      voxel_size = voxel_size * 3
    try:
      if len(voxel_size) != 3:
        raise ValueError
      [float(v) for v in voxel_size]  #pylint: disable=expression-not-assigned
    except ValueError as exception:
      raise MRtrixError('voxel size needs to be a single or three comma-separated floating point numbers; received: ' + str(app.ARGS.voxel_size)) from exception

  agg_measure = 'mean'
  if app.ARGS.aggregate is not None:
    if not app.ARGS.aggregate in AGGREGATION_MODES:
      app.error("aggregation type must be one of %s. provided: %s" % (str(AGGREGATION_MODES), app.ARGS.aggregate))
    agg_measure = app.ARGS.aggregate

  agg_weights = app.ARGS.aggregation_weights
  if agg_weights is not None:
    agg_measure = "weighted_" + agg_measure
    if agg_measure != 'weighted_mean':
      app.error("aggregation weights require '-aggregate mean' option. provided: %s" % (app.ARGS.aggregate))
      if not os.path.isfile(app.ARGS.aggregation_weights):
        app.error("aggregation weights file not found: %s" % app.ARGS.aggregation_weights)

  initial_alignment = app.ARGS.initial_alignment
  if initial_alignment not in ["mass", "robust_mass", "geometric", "none"]:
    raise MRtrixError('initial_alignment must be one of ' + " ".join(["mass", "robust_mass", "geometric", "none"]) + " provided: " + str(initial_alignment))

  linear_estimator = app.ARGS.linear_estimator
  if linear_estimator and not linear_estimator.lower() == 'none':
    if not dolinear:
      raise MRtrixError('linear_estimator specified when no linear registration is requested')
    if linear_estimator not in ["l1", "l2", "lp"]:
      raise MRtrixError('linear_estimator must be one of ' + " ".join(["l1", "l2", "lp"]) + " provided: " + str(linear_estimator))

  use_masks = False
  mask_files = []
  if app.ARGS.mask_dir:
    use_masks = True
    app.ARGS.mask_dir = relpath(app.ARGS.mask_dir)
    if not os.path.isdir(app.ARGS.mask_dir):
      raise MRtrixError('mask directory not found')
    mask_files = sorted(path.all_in_dir(app.ARGS.mask_dir, dir_path=False))
    if len(mask_files) < len(in_files[0]):
      raise MRtrixError('there are not enough mask images for the number of images in the input directory')

  if not use_masks:
    app.warn('no masks input. Use input masks to reduce computation time and improve robustness')

  if app.ARGS.template_mask and not use_masks:
    raise MRtrixError('you cannot output a template mask because no subject masks were input using -mask_dir')

  nanmask_input = app.ARGS.nanmask
  if nanmask_input and not use_masks:
    raise MRtrixError('you cannot use NaN masking when no subject masks were input using -mask_dir')

  ins, xcontrast_xsubject_pre_postfix = parse_input_files(in_files, mask_files, cns, agg_weights)

  leave_one_out = 'auto'
  if app.ARGS.leave_one_out is not None:
    leave_one_out = app.ARGS.leave_one_out
    if not leave_one_out in ['0', '1', 'auto']:
      raise MRtrixError('leave_one_out not understood: ' + str(leave_one_out))
  if leave_one_out == 'auto':
    leave_one_out = 2 < len(ins) < 15
  else:
    leave_one_out = bool(int(leave_one_out))
  if leave_one_out:
    app.console('performing leave-one-out registration')
    # check that at sum of weights is positive for any grouping if weighted aggregation is used
    weights = [float(inp.aggregation_weight) for inp in ins if inp.aggregation_weight is not None]
    if weights and sum(weights) - max(weights) <= 0:
      raise MRtrixError('leave-one-out registration requires positive aggregation weights in all groupings')

  noreorientation = app.ARGS.noreorientation

  do_pause_on_warn = True
  if app.ARGS.linear_no_pause:
    do_pause_on_warn = False
    if not dolinear:
      raise MRtrixError("linear option set when no linear registration is performed")

  if len(app.ARGS.template) != n_contrasts:
    raise MRtrixError('mismatch between number of output templates (%i) ' % len(app.ARGS.template) +
                      'and number of contrasts (%i)' % n_contrasts)
  for templ in app.ARGS.template:
    app.check_output_path(templ)

  if app.ARGS.warp_dir:
    app.ARGS.warp_dir = relpath(app.ARGS.warp_dir)
    app.check_output_path(app.ARGS.warp_dir)

  if app.ARGS.transformed_dir:
    app.ARGS.transformed_dir = [relpath(d) for d in app.ARGS.transformed_dir.split(',')]
    if len(app.ARGS.transformed_dir) != n_contrasts:
      raise MRtrixError('require multiple comma separated transformed directories if multi-contrast registration is used')
    for tdir in app.ARGS.transformed_dir:
      app.check_output_path(tdir)

  if app.ARGS.linear_transformations_dir:
    if not dolinear:
      raise MRtrixError("linear option set when no linear registration is performed")
    app.ARGS.linear_transformations_dir = relpath(app.ARGS.linear_transformations_dir)
    app.check_output_path(app.ARGS.linear_transformations_dir)

  # automatically detect SH series in each contrast
  do_fod_registration = False  # in any contrast
  cns.n_volumes = []
  cns.fod_reorientation = []
  for cid in range(n_contrasts):
    header = image.Header(ins[0].get_ims_path(False)[cid])
    image_size = header.size()
    if len(image_size) < 3 or len(image_size) > 4:
      raise MRtrixError('only 3 and 4 dimensional images can be used to build a template')
    if len(image_size) == 4:
      cns.fod_reorientation.append(header.is_sh() and not noreorientation)
      cns.n_volumes.append(image_size[3])
      do_fod_registration = do_fod_registration or cns.fod_reorientation[-1]
    else:
      cns.fod_reorientation.append(False)
      cns.n_volumes.append(0)
  if do_fod_registration:
    app.console("SH Series detected, performing FOD registration in contrast: " +
                ', '.join(app.ARGS.input_dir[cid] for cid in range(n_contrasts) if cns.fod_reorientation[cid]))
  c_mrtransform_reorientation = [' -reorient_fod ' + ('yes' if cns.fod_reorientation[cid] else 'no') + ' '
                                 for cid in range(n_contrasts)]

  if nanmask_input:
    app.console("NaN masking transformed images")

  # rigid options
  if app.ARGS.rigid_scale:
    rigid_scales = [float(x) for x in app.ARGS.rigid_scale.split(',')]
    if not dorigid:
      raise MRtrixError("rigid_scales option set when no rigid registration is performed")
  else:
    rigid_scales = DEFAULT_RIGID_SCALES
  if app.ARGS.rigid_lmax:
    if not dorigid:
      raise MRtrixError("rigid_lmax option set when no rigid registration is performed")
    rigid_lmax = [int(x) for x in app.ARGS.rigid_lmax.split(',')]
    if do_fod_registration and len(rigid_scales) != len(rigid_lmax):
      raise MRtrixError('rigid_scales and rigid_lmax schedules are not equal in length: scales stages: %s, lmax stages: %s' % (len(rigid_scales), len(rigid_lmax)))
  else:
    rigid_lmax = DEFAULT_RIGID_LMAX

  rigid_niter = [100] * len(rigid_scales)
  if app.ARGS.rigid_niter:
    if not dorigid:
      raise MRtrixError("rigid_niter specified when no rigid registration is performed")
    rigid_niter = [int(x) for x in app.ARGS.rigid_niter.split(',')]
    if len(rigid_niter) == 1:
      rigid_niter = rigid_niter * len(rigid_scales)
    elif len(rigid_scales) != len(rigid_niter):
      raise MRtrixError('rigid_scales and rigid_niter schedules are not equal in length: scales stages: %s, niter stages: %s' % (len(rigid_scales), len(rigid_niter)))

  # affine options
  if app.ARGS.affine_scale:
    affine_scales = [float(x) for x in app.ARGS.affine_scale.split(',')]
    if not doaffine:
      raise MRtrixError("affine_scale option set when no affine registration is performed")
  else:
    affine_scales = DEFAULT_AFFINE_SCALES
  if app.ARGS.affine_lmax:
    if not doaffine:
      raise MRtrixError("affine_lmax option set when no affine registration is performed")
    affine_lmax = [int(x) for x in app.ARGS.affine_lmax.split(',')]
    if do_fod_registration and len(affine_scales) != len(affine_lmax):
      raise MRtrixError('affine_scales and affine_lmax schedules are not equal in length: scales stages: %s, lmax stages: %s' % (len(affine_scales), len(affine_lmax)))
  else:
    affine_lmax = DEFAULT_AFFINE_LMAX

  affine_niter = [500] * len(affine_scales)
  if app.ARGS.affine_niter:
    if not doaffine:
      raise MRtrixError("affine_niter specified when no affine registration is performed")
    affine_niter = [int(x) for x in app.ARGS.affine_niter.split(',')]
    if len(affine_niter) == 1:
      affine_niter = affine_niter * len(affine_scales)
    elif len(affine_scales) != len(affine_niter):
      raise MRtrixError('affine_scales and affine_niter schedules are not equal in length: scales stages: %s, niter stages: %s' % (len(affine_scales), len(affine_niter)))

  linear_scales = []
  linear_lmax = []
  linear_niter = []
  linear_type = []
  if dorigid:
    linear_scales += rigid_scales
    linear_lmax += rigid_lmax
    linear_niter += rigid_niter
    linear_type += ['rigid'] * len(rigid_scales)

  if doaffine:
    linear_scales += affine_scales
    linear_lmax += affine_lmax
    linear_niter += affine_niter
    linear_type += ['affine'] * len(affine_scales)

  assert len(linear_type) == len(linear_scales)
  assert len(linear_scales) == len(linear_niter)
  if do_fod_registration:
    if len(linear_lmax) != len(linear_niter):
      mismatch = []
      if len(rigid_lmax) != len(rigid_niter):
        mismatch += ['rigid: lmax stages: %s, niter stages: %s' % (len(rigid_lmax), len(rigid_niter))]
      if len(affine_lmax) != len(affine_niter):
        mismatch += ['affine: lmax stages: %s, niter stages: %s' % (len(affine_lmax), len(affine_niter))]
      raise MRtrixError('linear registration: lmax and niter schedules are not equal in length: %s' % (', '.join(mismatch)))
  app.console('-' * 60)
  app.console('initial alignment of images: %s' % initial_alignment)
  app.console('-' * 60)
  if n_contrasts > 1:
    for cid in range(n_contrasts):
      app.console('\tcontrast "%s": %s, ' % (cns.suff[cid], cns.names[cid]) +
                  'objective weight: %s' % cns.mc_weight_initial_alignment[cid])

  if dolinear:
    app.console('-' * 60)
    app.console('linear registration stages:')
    app.console('-' * 60)
    if n_contrasts > 1:
      for cid in range(n_contrasts):
        msg = '\tcontrast "%s": %s' % (cns.suff[cid], cns.names[cid])
        if 'rigid' in linear_type:
          msg += ', objective weight rigid: %s' % cns.mc_weight_rigid[cid]
        if 'affine' in linear_type:
          msg += ', objective weight affine: %s' % cns.mc_weight_affine[cid]
        app.console(msg)

    if do_fod_registration:
      for istage, [tpe, scale, lmax, niter] in enumerate(zip(linear_type, linear_scales, linear_lmax, linear_niter)):
        app.console('(%02i) %s scale: %.4f, niter: %i, lmax: %i' % (istage, tpe.ljust(9), scale, niter, lmax))
    else:
      for istage, [tpe, scale, niter] in enumerate(zip(linear_type, linear_scales, linear_niter)):
        app.console('(%02i) %s scale: %.4f, niter: %i, no reorientation' % (istage, tpe.ljust(9), scale, niter))

  datatype_option = ' -datatype float32'
  outofbounds_option = ' -nan'

  if not dononlinear:
    nl_scales = []
    nl_lmax = []
    nl_niter = []
    if app.ARGS.warp_dir:
      raise MRtrixError('warp_dir specified when no nonlinear registration is performed')
  else:
    nl_scales = [float(x) for x in app.ARGS.nl_scale.split(',')] if app.ARGS.nl_scale else DEFAULT_NL_SCALES
    nl_niter = [int(x) for x in app.ARGS.nl_niter.split(',')] if app.ARGS.nl_niter else DEFAULT_NL_NITER
    nl_lmax = [int(x) for x in app.ARGS.nl_lmax.split(',')] if app.ARGS.nl_lmax else DEFAULT_NL_LMAX

    if len(nl_scales) != len(nl_niter):
      raise MRtrixError('nl_scales and nl_niter schedules are not equal in length: scales stages: %s, niter stages: %s' % (len(nl_scales), len(nl_niter)))

    app.console('-' * 60)
    app.console('nonlinear registration stages:')
    app.console('-' * 60)
    if n_contrasts > 1:
      for cid in range(n_contrasts):
        app.console('\tcontrast "%s": %s, objective weight: %s' % (cns.suff[cid], cns.names[cid], cns.mc_weight_nl[cid]))

    if do_fod_registration:
      if len(nl_scales) != len(nl_lmax):
        raise MRtrixError('nl_scales and nl_lmax schedules are not equal in length: scales stages: %s, lmax stages: %s' % (len(nl_scales), len(nl_lmax)))

    if do_fod_registration:
      for istage, [scale, lmax, niter] in enumerate(zip(nl_scales, nl_lmax, nl_niter)):
        app.console('(%02i) nonlinear scale: %.4f, niter: %i, lmax: %i' % (istage, scale, niter, lmax))
    else:
      for istage, [scale, niter] in enumerate(zip(nl_scales, nl_niter)):
        app.console('(%02i) nonlinear scale: %.4f, niter: %i, no reorientation' % (istage, scale, niter))

  app.console('-' * 60)
  app.console('input images:')
  app.console('-' * 60)
  for inp in ins:
    app.console('\t' + inp.info())

  app.make_scratch_dir()
  app.goto_scratch_dir()

  for contrast in cns.suff:
    path.make_dir('input_transformed' + contrast)

  for contrast in cns.suff:
    path.make_dir('isfinite' + contrast)

  path.make_dir('linear_transforms_initial')
  path.make_dir('linear_transforms')
  for level in range(0, len(linear_scales)):
    path.make_dir('linear_transforms_%02i' % level)
  for level in range(0, len(nl_scales)):
    path.make_dir('warps_%02i' % level)

  if use_masks:
    path.make_dir('mask_transformed')
  write_log = (app.VERBOSITY >= 2)
  if write_log:
    path.make_dir('log')

  if initial_alignment == 'robust_mass':
    if not use_masks:
      raise MRtrixError('robust_mass initial alignment requires masks')
    path.make_dir('robust')

  if app.ARGS.copy_input:
    app.console('Copying images into scratch directory')
    for inp in ins:
      inp.cache_local()

  # Make initial template in average space using first contrast
  app.console('Generating initial template')
  input_filenames = [inp.get_ims_path(False)[0] for inp in ins]
  if voxel_size is None:
    run.command(['mraverageheader', input_filenames, 'average_header.mif', '-fill'])
  else:
    run.command(['mraverageheader', '-fill', input_filenames, '-', '|',
                 'mrgrid', '-', 'regrid', '-voxel', ','.join(map(str, voxel_size)), 'average_header.mif'])

  # crop average space to extent defined by original masks
  if use_masks:
    progress = app.ProgressBar('Importing input masks to average space for template cropping', len(ins))
    for inp in ins:
      run.command('mrtransform ' + inp.msk_path + ' -interp nearest -template average_header.mif ' + inp.msk_transformed)
      progress.increment()
    progress.done()
    run.command(['mrmath', [inp.msk_transformed for inp in ins], 'max', 'mask_initial.mif'])
    run.command('mrgrid average_header.mif crop -mask mask_initial.mif average_header_cropped.mif')
    run.function(os.remove, 'mask_initial.mif')
    run.function(os.remove, 'average_header.mif')
    run.function(shutil.move, 'average_header_cropped.mif', 'average_header.mif')
    progress = app.ProgressBar('Erasing temporary mask images', len(ins))
    for inp in ins:
      run.function(os.remove, inp.msk_transformed)
      progress.increment()
    progress.done()

  # create average space headers for other contrasts
  if n_contrasts > 1:
    avh3d = 'average_header3d.mif'
    avh4d = 'average_header4d.mif'
    if len(image.Header('average_header.mif').size()) == 3:
      run.command('mrconvert average_header.mif ' + avh3d)
    else:
      run.command('mrconvert average_header.mif -coord 3 0 -axes 0,1,2 ' + avh3d)
    run.command('mrconvert ' + avh3d + ' -axes 0,1,2,-1 ' + avh4d)
    for cid in range(n_contrasts):
      if cns.n_volumes[cid] == 0:
        run.function(copy, avh3d, 'average_header' + cns.suff[cid] + '.mif')
      elif cns.n_volumes[cid] == 1:
        run.function(copy, avh4d, 'average_header' + cns.suff[cid] + '.mif')
      else:
        run.command(['mrcat', [avh3d] * cns.n_volumes[cid], '-axis', '3', 'average_header' + cns.suff[cid] + '.mif'])
    run.function(os.remove, avh3d)
    run.function(os.remove, avh4d)
  else:
    run.function(shutil.move, 'average_header.mif', 'average_header' + cns.suff[0] + '.mif')

  cns.templates = ['average_header' + csuff + '.mif' for csuff in cns.suff]

  if initial_alignment == 'none':
    progress = app.ProgressBar('Resampling input images to template space with no initial alignment', len(ins) * n_contrasts)
    for inp in ins:
      for cid in range(n_contrasts):
        run.command('mrtransform ' + inp.ims_path[cid] + c_mrtransform_reorientation[cid] + ' -interp linear ' +
                    '-template ' + cns.templates[cid] + ' ' + inp.ims_transformed[cid] +
                    outofbounds_option +
                    datatype_option)
        progress.increment()
    progress.done()

    if use_masks:
      progress = app.ProgressBar('Reslicing input masks to average header', len(ins))
      for inp in ins:
        run.command('mrtransform ' + inp.msk_path + ' ' + inp.msk_transformed + ' ' +
                    '-interp nearest -template ' + cns.templates[0] + ' ' +
                    datatype_option)
        progress.increment()
      progress.done()

    if nanmask_input:
      inplace_nan_mask([inp.ims_transformed[cid] for inp in ins for cid in range(n_contrasts)],
                       [inp.msk_transformed for inp in ins for cid in range(n_contrasts)])

    if leave_one_out:
      calculate_isfinite(ins, cns)

    if not dolinear:
      for inp in ins:
        with open(os.path.join('linear_transforms_initial', inp.uid + '.txt'), 'w', encoding='utf-8') as fout:
          fout.write('1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n')

    run.function(copy, 'average_header' + cns.suff[0] + '.mif', 'average_header.mif')

  else:
    progress = app.ProgressBar('Performing initial rigid registration to template', len(ins))
    mask_option = ''
    cid = 0
    lmax_option = ' -rigid_lmax 0 ' if cns.fod_reorientation[cid] else ' -noreorientation '
    contrast_weight_option = cns.initial_alignment_weight_option
    for inp in ins:
      output_option = ' -rigid ' + os.path.join('linear_transforms_initial', inp.uid + '.txt')
      images = ' '.join([p + ' ' + t for p, t in zip(inp.ims_path, cns.templates)])
      if use_masks:
        mask_option = ' -mask1 ' + inp.msk_path
        if initial_alignment == 'robust_mass':
          if not os.path.isfile('robust/template.mif'):
            if cns.n_volumes[cid] > 0:
              run.command('mrconvert ' + cns.templates[cid] + ' -coord 3 0 - | mrconvert - -axes 0,1,2 robust/template.mif')
            else:
              run.command('mrconvert ' + cns.templates[cid] + ' robust/template.mif')
          if n_contrasts > 1:
            cmd = ['mrcalc', inp.ims_path[cid], cns.mc_weight_initial_alignment[cid], '-mult']
            for cid in range(1, n_contrasts):
              cmd += [inp.ims_path[cid], cns.mc_weight_initial_alignment[cid], '-mult', '-add']
            contrast_weight_option = ''
            run.command(' '.join(cmd) +
                        ' - | mrfilter - zclean -zlower 3 -zupper 3 robust/image_' + inp.uid + '.mif'
                        ' -maskin ' + inp.msk_path + ' -maskout robust/mask_' + inp.uid + '.mif')
          else:
            run.command('mrfilter ' + inp.ims_path[0] + ' zclean -zlower 3 -zupper 3 robust/image_' + inp.uid + '.mif' +
                        ' -maskin ' + inp.msk_path + ' -maskout robust/mask_' + inp.uid + '.mif')
          images = 'robust/image_' + inp.uid + '.mif robust/template.mif'
          mask_option = ' -mask1 ' + 'robust/mask_' + inp.uid + '.mif'
          lmax_option = ''

      run.command('mrregister ' + images +
                  mask_option +
                  ' -rigid_scale 1 ' +
                  ' -rigid_niter 0 ' +
                  ' -type rigid ' +
                  lmax_option +
                  contrast_weight_option +
                  ' -rigid_init_translation ' + initial_alignment.replace('robust_', '') + ' ' +
                  datatype_option +
                  output_option)
      # translate input images to centre of mass without interpolation
      for cid in range(n_contrasts):
        run.command('mrtransform ' + inp.ims_path[cid] + c_mrtransform_reorientation[cid] +
                    ' -linear ' + os.path.join('linear_transforms_initial', inp.uid + '.txt') +
                    ' ' + inp.ims_transformed[cid] + "_translated.mif" + datatype_option)
      if use_masks:
        run.command('mrtransform ' + inp.msk_path +
                    ' -linear ' + os.path.join('linear_transforms_initial', inp.uid + '.txt') +
                    ' ' + inp.msk_transformed + "_translated.mif" +
                    datatype_option)
      progress.increment()
    # update average space of first contrast to new extent, delete other average space images
    run.command(['mraverageheader', [inp.ims_transformed[cid] + '_translated.mif' for inp in ins], 'average_header_tight.mif'])
    progress.done()

    if voxel_size is None:
      run.command('mrgrid average_header_tight.mif pad -uniform 10 average_header.mif', force=True)
    else:
      run.command('mrgrid average_header_tight.mif pad -uniform 10 - | '
                  'mrgrid - regrid -voxel ' + ','.join(map(str, voxel_size)) + ' average_header.mif', force=True)
    run.function(os.remove, 'average_header_tight.mif')
    for cid in range(1, n_contrasts):
      run.function(os.remove, 'average_header' + cns.suff[cid] + '.mif')

    if use_masks:
      # reslice masks
      progress = app.ProgressBar('Reslicing input masks to average header', len(ins))
      for inp in ins:
        run.command('mrtransform ' + inp.msk_transformed + '_translated.mif' + ' ' + inp.msk_transformed + ' ' +
                    '-interp nearest -template average_header.mif' + datatype_option)
        progress.increment()
      progress.done()
      # crop average space to extent defined by translated masks
      run.command(['mrmath', [inp.msk_transformed for inp in ins], 'max', 'mask_translated.mif'])
      run.command('mrgrid average_header.mif crop -mask mask_translated.mif average_header_cropped.mif')
      # pad average space to allow for deviation from initial alignment
      run.command('mrgrid average_header_cropped.mif pad -uniform 10 average_header.mif', force=True)
      run.function(os.remove, 'average_header_cropped.mif')
      # reslice masks
      progress = app.ProgressBar('Reslicing masks to new padded average header', len(ins))
      for inp in ins:
        run.command('mrtransform ' + inp.msk_transformed + '_translated.mif ' + inp.msk_transformed + ' ' +
                    '-interp nearest -template average_header.mif' + datatype_option, force=True)
        run.function(os.remove, inp.msk_transformed + '_translated.mif')
        progress.increment()
      progress.done()
      run.function(os.remove, 'mask_translated.mif')

    # reslice images
    progress = app.ProgressBar('Reslicing input images to average header', len(ins) * n_contrasts)
    for cid in range(n_contrasts):
      for inp in ins:
        run.command('mrtransform ' + c_mrtransform_reorientation[cid] + inp.ims_transformed[cid] + '_translated.mif ' +
                    inp.ims_transformed[cid] + ' ' +
                    ' -interp linear -template average_header.mif' +
                    outofbounds_option +
                    datatype_option)
        run.function(os.remove, inp.ims_transformed[cid] + '_translated.mif')
        progress.increment()
    progress.done()

    if nanmask_input:
      inplace_nan_mask([inp.ims_transformed[cid] for inp in ins for cid in range(n_contrasts)],
                       [inp.msk_transformed for inp in ins for cid in range(n_contrasts)])

    if leave_one_out:
      calculate_isfinite(ins, cns)

  cns.templates = ['initial_template' + contrast + '.mif' for contrast in cns.suff]
  for cid in range(n_contrasts):
    aggregate(ins, 'initial_template' + cns.suff[cid] + '.mif', cid, agg_measure)
    if cns.n_volumes[cid] == 1:
      run.function(shutil.move, 'initial_template' + cns.suff[cid] + '.mif', 'tmp.mif')
      run.command('mrconvert tmp.mif initial_template' + cns.suff[cid] + '.mif -axes 0,1,2,-1')

  # Optimise template with linear registration
  if not dolinear:
    for inp in ins:
      run.function(copy, os.path.join('linear_transforms_initial', inp.uid+'.txt'),
                   os.path.join('linear_transforms', inp.uid+'.txt'))
  else:
    level = 0
    regtype = linear_type[0]
    def linear_msg():
      return 'Optimising template with linear registration (stage {0} of {1}; {2})'.format(level + 1, len(linear_scales), regtype)
    progress = app.ProgressBar(linear_msg, len(linear_scales) * len(ins) * (1 + n_contrasts + int(use_masks)))
    for level, (regtype, scale, niter, lmax) in enumerate(zip(linear_type, linear_scales, linear_niter, linear_lmax)):
      for inp in ins:
        initialise_option = ''
        if use_masks:
          mask_option = ' -mask1 ' + inp.msk_path
        else:
          mask_option = ''
        lmax_option = ' -noreorientation'
        metric_option = ''
        mrregister_log_option = ''
        if regtype == 'rigid':
          scale_option = ' -rigid_scale ' + str(scale)
          niter_option = ' -rigid_niter ' + str(niter)
          regtype_option = ' -type rigid'
          output_option = ' -rigid ' + os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt')
          contrast_weight_option = cns.rigid_weight_option
          initialise_option = (' -rigid_init_matrix ' +
                               os.path.join('linear_transforms_%02i' % (level - 1) if level > 0 else 'linear_transforms_initial', inp.uid + '.txt'))
          if do_fod_registration:
            lmax_option = ' -rigid_lmax ' + str(lmax)
          if linear_estimator:
            metric_option = ' -rigid_metric.diff.estimator ' + linear_estimator
          if app.VERBOSITY >= 2:
            mrregister_log_option = ' -info -rigid_log ' + os.path.join('log', inp.uid + contrast[cid] + "_" + str(level) + '.log')
        else:
          scale_option = ' -affine_scale ' + str(scale)
          niter_option = ' -affine_niter ' + str(niter)
          regtype_option = ' -type affine'
          output_option = ' -affine ' + os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt')
          contrast_weight_option = cns.affine_weight_option
          initialise_option = (' -affine_init_matrix ' +
                               os.path.join('linear_transforms_%02i' % (level - 1) if level > 0 else 'linear_transforms_initial', inp.uid + '.txt'))
          if do_fod_registration:
            lmax_option = ' -affine_lmax ' + str(lmax)
          if linear_estimator:
            metric_option = ' -affine_metric.diff.estimator ' + linear_estimator
          if write_log:
            mrregister_log_option = ' -info -affine_log ' + os.path.join('log', inp.uid + contrast[cid] + "_" + str(level) + '.log')

        if leave_one_out:
          tmpl = []
          for cid in range(n_contrasts):
            isfinite = 'isfinite%s/%s.mif' % (cns.suff[cid], inp.uid)
            weight = inp.aggregation_weight if inp.aggregation_weight is not None else '1'
            # loo = (template * weighted sum - weight * this) / (weighted sum - weight)
            run.command('mrcalc ' + cns.isfinite_count[cid] + ' ' + isfinite + ' -sub - | mrcalc ' + cns.templates[cid] +
                        ' ' + cns.isfinite_count[cid] + ' -mult ' + inp.ims_transformed[cid] + ' ' + weight + ' -mult ' +
                        ' -sub - -div loo_%s' % cns.templates[cid], force=True)
            tmpl.append('loo_%s' % cns.templates[cid])
          images = ' '.join([p + ' ' + t for p, t in zip(inp.ims_path, tmpl)])
        else:
          images = ' '.join([p + ' ' + t for p, t in zip(inp.ims_path, cns.templates)])
        command = 'mrregister ' + images + \
                  initialise_option + \
                  mask_option + \
                  scale_option + \
                  niter_option + \
                  lmax_option + \
                  regtype_option + \
                  metric_option + \
                  datatype_option + \
                  contrast_weight_option + \
                  output_option + \
                  mrregister_log_option
        run.command(command, force=True)
        check_linear_transformation(os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt'), command,
                                    pause_on_warn=do_pause_on_warn)
        if leave_one_out:
          for im_temp in tmpl:
            run.function(os.remove, im_temp)
        progress.increment()

      # Here we ensure the overall template properties don't change (too much) over levels
      # The reference is the initialisation as that's used to construct the average space.
      # T_i: linear trafo for case i, i.e. template(x) = E [ image_i(T_i x) ]
      # R_i: inital trafo for case i (identity if initial alignment is none)
      # A = E[ T_i ]: average of current trafos
      # B = E[ R_i ]: average of initial trafos
      # C_i' = T_i B A^{-1}: "drift" corrected T_i
      # T_i <- C_i
      # Notes:
      #   - This approximately stabilises E[ T_i ], its' relatively close to B
      #   - Not sure whether it's preferable to stabilise E[ T_i^{-1} ]
      #   - If one subject's registration fails, this will affect the average and therefore the template which could result in instable behaviour.
      #   - The template appearance changes slightly over levels, but the template and trafos are affected in the same way so should not affect template convergence.
      if not app.ARGS.linear_no_drift_correction:
        run.command(['transformcalc', [os.path.join('linear_transforms_initial', inp.uid + '.txt') for _inp in ins],
                    'average', 'linear_transform_average_init.txt', '-quiet'], force=True)
        run.command(['transformcalc', [os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt') for _inp in ins],
                    'average', 'linear_transform_average_%02i_uncorrected.txt' % level, '-quiet'], force=True)
        run.command(['transformcalc', 'linear_transform_average_%02i_uncorrected.txt' % level,
                    'invert', 'linear_transform_average_%02i_uncorrected_inv.txt' % level, '-quiet'], force=True)

        transform_average_init = matrix.load_transform('linear_transform_average_init.txt')
        transform_average_current_inv = matrix.load_transform('linear_transform_average_%02i_uncorrected_inv.txt' % level)

        transform_update = matrix.dot(transform_average_init, transform_average_current_inv)
        matrix.save_transform(os.path.join('linear_transforms_%02i_drift_correction.txt' %  level), transform_update, force=True)
        if regtype == 'rigid':
          run.command('transformcalc ' + os.path.join('linear_transforms_%02i_drift_correction.txt' %  level) +
                      ' rigid ' + os.path.join('linear_transforms_%02i_drift_correction.txt' %  level) + ' -quiet', force=True)
          transform_update = matrix.load_transform(os.path.join('linear_transforms_%02i_drift_correction.txt' %  level))

        for inp in ins:
          transform = matrix.load_transform('linear_transforms_%02i/' % level + inp.uid + '.txt')
          transform_updated = matrix.dot(transform, transform_update)
          run.function(copy, 'linear_transforms_%02i/' % level + inp.uid + '.txt', 'linear_transforms_%02i/' % level + inp.uid + '.precorrection')
          matrix.save_transform(os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt'), transform_updated, force=True)

        # compute average trafos and its properties for easier debugging
        run.command(['transformcalc', [os.path.join('linear_transforms_%02i' % level, _inp.uid + '.txt') for _inp in ins],
                    'average', 'linear_transform_average_%02i.txt' % level, '-quiet'], force=True)
        run.command('transformcalc linear_transform_average_%02i.txt decompose linear_transform_average_%02i.dec' % (level, level), force=True)


      for cid in range(n_contrasts):
        for inp in ins:
          run.command('mrtransform ' + c_mrtransform_reorientation[cid] + inp.ims_path[cid] +
                      ' -template ' + cns.templates[cid] +
                      ' -linear ' + os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt') +
                      ' ' + inp.ims_transformed[cid] +
                      outofbounds_option +
                      datatype_option,
                      force=True)
          progress.increment()

      if use_masks:
        for inp in ins:
          run.command('mrtransform ' + inp.msk_path +
                      ' -template ' + cns.templates[0] +
                      ' -interp nearest' +
                      ' -linear ' + os.path.join('linear_transforms_%02i' % level, inp.uid + '.txt') +
                      ' ' + inp.msk_transformed,
                      force=True)
          progress.increment()

      if nanmask_input:
        inplace_nan_mask([inp.ims_transformed[cid] for inp in ins for cid in range(n_contrasts)],
                         [inp.msk_transformed for inp in ins for cid in range(n_contrasts)])

      if leave_one_out:
        calculate_isfinite(ins, cns)

      for cid in range(n_contrasts):
        if level > 0 and app.ARGS.delete_temporary_files:
          os.remove(cns.templates[cid])
        cns.templates[cid] = 'linear_template%02i%s.mif' % (level, cns.suff[cid])
        aggregate(ins, cns.templates[cid], cid, agg_measure)
        if cns.n_volumes[cid] == 1:
          run.function(shutil.move, cns.templates[cid], 'tmp.mif')
          run.command('mrconvert tmp.mif ' + cns.templates[cid] + ' -axes 0,1,2,-1')
          run.function(os.remove, 'tmp.mif')

    for entry in os.listdir('linear_transforms_%02i' % level):
      run.function(copy, os.path.join('linear_transforms_%02i' % level, entry), os.path.join('linear_transforms', entry))
    progress.done()

  # Create a template mask for nl registration by taking the intersection of all transformed input masks and dilating
  if use_masks and (dononlinear or app.ARGS.template_mask):
    run.command(['mrmath', path.all_in_dir('mask_transformed')] +
                'min - | maskfilter - median - | maskfilter - dilate -npass 5 init_nl_template_mask.mif'.split(), force=True)
    current_template_mask = 'init_nl_template_mask.mif'

  if dononlinear:
    path.make_dir('warps')
    level = 0
    def nonlinear_msg():
      return 'Optimising template with non-linear registration (stage {0} of {1})'.format(level + 1, len(nl_scales))
    progress = app.ProgressBar(nonlinear_msg, len(nl_scales) * len(ins))
    for level, (scale, niter, lmax) in enumerate(zip(nl_scales, nl_niter, nl_lmax)):
      for inp in ins:
        if level > 0:
          initialise_option = ' -nl_init ' + os.path.join('warps_%02i' % (level - 1), inp.uid + '.mif')
          scale_option = ''
        else:
          scale_option = ' -nl_scale ' + str(scale)
          if not doaffine:  # rigid or no previous linear stage
            initialise_option = ' -rigid_init_matrix ' + os.path.join('linear_transforms', inp.uid + '.txt')
          else:
            initialise_option = ' -affine_init_matrix ' + os.path.join('linear_transforms', inp.uid + '.txt')

        if use_masks:
          mask_option = ' -mask1 ' + inp.msk_path + ' -mask2 ' + current_template_mask
        else:
          mask_option = ''

        if do_fod_registration:
          lmax_option = ' -nl_lmax ' + str(lmax)
        else:
          lmax_option = ' -noreorientation'

        contrast_weight_option = cns.nl_weight_option

        if leave_one_out:
          tmpl = []
          for cid in range(n_contrasts):
            isfinite = 'isfinite%s/%s.mif' % (cns.suff[cid], inp.uid)
            weight = inp.aggregation_weight if inp.aggregation_weight is not None else '1'
            # loo = (template * weighted sum - weight * this) / (weighted sum - weight)
            run.command('mrcalc ' + cns.isfinite_count[cid] + ' ' + isfinite + ' -sub - | mrcalc ' + cns.templates[cid] +
                        ' ' + cns.isfinite_count[cid] + ' -mult ' + inp.ims_transformed[cid] + ' ' + weight + ' -mult ' +
                        ' -sub - -div loo_%s' % cns.templates[cid], force=True)
            tmpl.append('loo_%s' % cns.templates[cid])
          images = ' '.join([p + ' ' + t for p, t in zip(inp.ims_path, tmpl)])
        else:
          images = ' '.join([p + ' ' + t for p, t in zip(inp.ims_path, cns.templates)])
        run.command('mrregister ' + images +
                    ' -type nonlinear' +
                    ' -nl_niter ' + str(nl_niter[level]) +
                    ' -nl_warp_full ' + os.path.join('warps_%02i' % level, inp.uid + '.mif') +
                    ' -transformed ' +
                    ' -transformed '.join([inp.ims_transformed[cid] for cid in range(n_contrasts)]) + ' ' +
                    ' -nl_update_smooth ' + app.ARGS.nl_update_smooth +
                    ' -nl_disp_smooth ' + app.ARGS.nl_disp_smooth +
                    ' -nl_grad_step ' + app.ARGS.nl_grad_step +
                    initialise_option +
                    contrast_weight_option +
                    scale_option +
                    mask_option +
                    datatype_option +
                    outofbounds_option +
                    lmax_option,
                    force=True)

        if use_masks:
          run.command('mrtransform ' + inp.msk_path +
                      ' -template ' + cns.templates[0] +
                      ' -warp_full ' + os.path.join('warps_%02i' % level, inp.uid + '.mif') +
                      ' ' + inp.msk_transformed +
                      ' -interp nearest ',
                      force=True)

        if leave_one_out:
          for im_temp in tmpl:
            run.function(os.remove, im_temp)

        if level > 0:
          run.function(os.remove, os.path.join('warps_%02i' % (level - 1), inp.uid + '.mif'))

        progress.increment(nonlinear_msg())

      if nanmask_input:
        inplace_nan_mask([_inp.ims_transformed[cid] for _inp in ins for cid in range(n_contrasts)],
                         [_inp.msk_transformed for _inp in ins for cid in range(n_contrasts)])

      if leave_one_out:
        calculate_isfinite(ins, cns)

      for cid in range(n_contrasts):
        if level > 0 and app.ARGS.delete_temporary_files:
          os.remove(cns.templates[cid])
        cns.templates[cid] = 'nl_template%02i%s.mif' % (level, cns.suff[cid])
        aggregate(ins, cns.templates[cid], cid, agg_measure)
        if cns.n_volumes[cid] == 1:
          run.function(shutil.move, cns.templates[cid], 'tmp.mif')
          run.command('mrconvert tmp.mif ' + cns.templates[cid] + ' -axes 0,1,2,-1')
          run.function(os.remove, 'tmp.mif')

      if use_masks:
        run.command(['mrmath', path.all_in_dir('mask_transformed')] +
                    'min - | maskfilter - median - | '.split() +
                    ('maskfilter - dilate -npass 5 nl_template_mask' + str(level) + '.mif').split())
        current_template_mask = 'nl_template_mask' + str(level) + '.mif'

      if level < len(nl_scales) - 1:
        if scale < nl_scales[level + 1]:
          upsample_factor = nl_scales[level + 1] / scale
          for inp in ins:
            run.command('mrgrid ' + os.path.join('warps_%02i' % level, inp.uid + '.mif') +
                        ' regrid -scale %f tmp.mif' % upsample_factor, force=True)
            run.function(shutil.move, 'tmp.mif', os.path.join('warps_%02i' % level, inp.uid + '.mif'))
      else:
        for inp in ins:
          run.function(shutil.move, os.path.join('warps_%02i' % level, inp.uid + '.mif'), 'warps')
    progress.done()

  for cid in range(n_contrasts):
    run.command('mrconvert ' + cns.templates[cid] + ' ' + cns.templates_out[cid],
                mrconvert_keyval='NULL', force=app.FORCE_OVERWRITE)

  if app.ARGS.warp_dir:
    warp_path = path.from_user(app.ARGS.warp_dir, False)
    if os.path.exists(warp_path):
      run.function(shutil.rmtree, warp_path)
    os.makedirs(warp_path)
    progress = app.ProgressBar('Copying non-linear warps to output directory "' + warp_path + '"', len(ins))
    for inp in ins:
      keyval = image.Header(os.path.join('warps', inp.uid + '.mif')).keyval()
      keyval = dict((k, keyval[k]) for k in ('linear1', 'linear2'))
      json_path = os.path.join('warps', inp.uid + '.json')
      with open(json_path, 'w', encoding='utf-8') as json_file:
        json.dump(keyval, json_file)
      run.command('mrconvert ' + os.path.join('warps', inp.uid + '.mif') + ' ' +
                  shlex.quote(os.path.join(warp_path, xcontrast_xsubject_pre_postfix[0] +
                                           inp.uid + xcontrast_xsubject_pre_postfix[1] + '.mif')),
                  mrconvert_keyval=json_path, force=app.FORCE_OVERWRITE)
      progress.increment()
    progress.done()

  if app.ARGS.linear_transformations_dir:
    linear_transformations_path = path.from_user(app.ARGS.linear_transformations_dir, False)
    if os.path.exists(linear_transformations_path):
      run.function(shutil.rmtree, linear_transformations_path)
    os.makedirs(linear_transformations_path)
    for inp in ins:
      trafo = matrix.load_transform(os.path.join('linear_transforms', inp.uid + '.txt'))
      matrix.save_transform(os.path.join(linear_transformations_path,
                                         xcontrast_xsubject_pre_postfix[0] + inp.uid
                                         + xcontrast_xsubject_pre_postfix[1] + '.txt'),
                            trafo,
                            force=app.FORCE_OVERWRITE)

  if app.ARGS.transformed_dir:
    for cid, trdir in enumerate(app.ARGS.transformed_dir):
      transformed_path = path.from_user(trdir, False)
      if os.path.exists(transformed_path):
        run.function(shutil.rmtree, transformed_path)
      os.makedirs(transformed_path)
      progress = app.ProgressBar('Copying transformed images to output directory "' + transformed_path + '"', len(ins))
      for inp in ins:
        run.command(['mrconvert', inp.ims_transformed[cid], os.path.join(transformed_path, inp.ims_filenames[cid])],
                    mrconvert_keyval=inp.get_ims_path(False)[cid], force=app.FORCE_OVERWRITE)
        progress.increment()
      progress.done()

  if app.ARGS.template_mask:
    run.command('mrconvert ' + current_template_mask + ' ' + path.from_user(app.ARGS.template_mask, True),
                mrconvert_keyval='NULL', force=app.FORCE_OVERWRITE)

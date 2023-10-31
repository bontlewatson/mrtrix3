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

def usage(cmdline): #pylint: disable=unused-variable
  cmdline.set_author('Remika Mito (remika.mito@florey.edu.au) and Robert E. Smith (robert.smith@florey.edu.au)')
  cmdline.set_synopsis('Create a glass brain from mask input')
  cmdline.add_description('The output of this command is a glass brain image, which can be viewed '
                          'using the volume render option in mrview, and used for visualisation purposes to view results in 3D.')
  cmdline.add_description('While the name of this script indicates that a binary mask image is required as input, it can '
                          'also operate on a floating-point image. One way in which this can be exploited is to compute the mean '
                          'of all subject masks within template space, in which case this script will produce a smoother result '
                          'than if a binary template mask were to be used as input.')
  cmdline.add_argument('input', help='The input mask image')
  cmdline.add_argument('output', help='The output glass brain image')
  cmdline.add_argument('-dilate', type=int, default=2, help='Provide number of passes for dilation step; default = 2')
  cmdline.add_argument('-scale', type=float, default=2.0, help='Provide resolution upscaling value; default = 2.0')
  cmdline.add_argument('-smooth', type=float, default=1.0, help='Provide standard deviation of smoothing (in mm); default = 1.0')


def execute(): #pylint: disable=unused-variable
  from mrtrix3 import app, image, path, run #pylint: disable=no-name-in-module, import-outside-toplevel

  app.check_output_path(app.ARGS.output)

  # import data to scratch directory
  app.make_scratch_dir()
  run.command('mrconvert ' + path.from_user(app.ARGS.input) + ' ' + path.to_scratch('in.mif'))
  app.goto_scratch_dir()

  dilate_option = ' -npass ' + str(app.ARGS.dilate)
  scale_option = ' -scale ' + str(app.ARGS.scale)
  smooth_option = ' -stdev ' + str(app.ARGS.smooth)
  threshold_option = ' -abs 0.5'

  # check whether threshold should be fixed at 0.5 or computed automatically from the data
  if image.Header('in.mif').datatype() == 'Bit':
    app.debug('Input image is bitwise; no need to check image intensities')
  else:
    app.debug('Input image is not bitwise; checking distribution of image intensities')
    result_stat = image.statistics('in.mif')
    if not (result_stat.min == 0.0 and result_stat.max == 1.0):
      app.warn('Input image contains values outside of range [0.0, 1.0]; threshold will not be 0.5, but will instead be determined from the image data')
      threshold_option = ''
    else:
      app.debug('Input image values reside within [0.0, 1.0] range; fixed threshold of 0.5 will be used')

  # run upscaling step
  run.command('mrgrid in.mif regrid upsampled.mif' + scale_option)

  # run smoothing step
  run.command('mrfilter upsampled.mif smooth upsampled_smooth.mif' + smooth_option)

  # threshold image
  run.command('mrthreshold upsampled_smooth.mif upsampled_smooth_thresh.mif' + threshold_option)

  # dilate image for subtraction
  run.command('maskfilter upsampled_smooth_thresh.mif dilate upsampled_smooth_thresh_dilate.mif' + dilate_option)

  # create border
  run.command('mrcalc upsampled_smooth_thresh_dilate.mif upsampled_smooth_thresh.mif -xor out.mif -datatype bit')

  # create output image
  run.command('mrconvert out.mif ' + path.from_user(app.ARGS.output),
              mrconvert_keyval=path.from_user(app.ARGS.input, False),
              force=app.FORCE_OVERWRITE)

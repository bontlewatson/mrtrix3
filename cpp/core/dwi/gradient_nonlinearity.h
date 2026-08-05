/* Copyright (c) 2008-2026 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Covered Software is provided under this License on an "as is"
 * basis, without warranty of any kind, either expressed, implied, or
 * statutory, including, without limitation, warranties that the
 * Covered Software is free of defects, merchantable, fit for a
 * particular purpose or non-infringing.
 * See the Mozilla Public License v. 2.0 for more details.
 *
 * For more details, see http://www.mrtrix.org/.
 */

#pragma once

#include <Eigen/Dense>

#include "image.h"
#include "adapter/reslice.h"
#include "interp/linear.h"

namespace MR::DWI {

  class GradientNonLinearityCorrection {
    public:
    GradientNonLinearityCorrection(const Image<float>& grad_dev, const Header& dwi_header) : 
      reslicer (grad_dev, dwi_header) {}

      // compute the L tensor from grad_dev image at a voxel pos
      void compute_L(Eigen::Matrix3d &L, const Eigen::Vector3i &vox) {
        // assuming L(x) stored as [Lxx, Lxy, Lxz, Lyx, Lyy, Lyz, Lzx, Lzy, Lzz]

        reslicer.index(0) = vox[0];
        reslicer.index(1) = vox[1];
        reslicer.index(2) = vox[2];

        reslicer.index(3) = 0;
        for (int i = 0; i < 3; ++i) {
          for (int j =0; j < 3; ++j) {
            ++reslicer.index(3);
            L(i,j) = static_cast<double>(reslicer.value());
          }
        }
      }

      // compute the corrected gradient information & bvalues at a voxel pos
      void correct_grad(const Eigen::MatrixXd& grad, Eigen::MatrixXd &grad_corr, const Eigen::Matrix3d &L) {
        assert(grad_corr.size() == grad.size());
        // transformation
        const Eigen::Matrix3d IL = Eigen::Matrix3d::Identity()+L;
        const Eigen::Matrix3d affine = IL*flipMat();

        for (int N = 0; N < grad.rows(); ++N) {
          const double b = grad(N,3);
          const Eigen::Vector3d bv = grad.row(N).head<3>();
          // account for left-handed system, i.e. flip the x-axis
          const Eigen::Vector3d g = affine * bv;
          const double g_norm = g.norm();

          if (g_norm > 0.0) {
            grad_corr.row(N).head<3>() = g / g_norm;
            grad_corr(N,3) = b * g_norm * g_norm;
          } else {
            grad_corr.row(N).setZero();
          }
        }
      }

      // flip matrix for -ve x-axis:
      static const Eigen::Matrix3d& flipMat() {
        static const Eigen::Matrix3d mat = (Eigen::Matrix3d() <<
            1, -1, -1,
            -1,  1,  1,
            -1,  1,  1).finished();
        return mat;
      }

    private:
      Adapter::Reslice<Interp::Linear, Image<float> > reslicer;

  };

}


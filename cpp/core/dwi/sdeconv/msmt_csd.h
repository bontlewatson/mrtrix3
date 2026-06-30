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
#include <optional>

#include "eigen_plugins/eigen_plugins.h"
#include "dwi/directions/predefined.h"
#include "dwi/directions/validate.h"
#include "dwi/gradient.h"
#include "dwi/sdeconv/se_response.h"
#include "dwi/shells.h"
#include "file/matrix.h"
#include "header.h"
#include "image.h"
#include "math/SH.h"
#include "math/ZSH.h"
#include "math/constrained_least_squares.h"
#include "math/math.h"
#include "math/sphere.h"
#include "types.h"

namespace MR::DWI::SDeconv {

constexpr uint32_t default_msmt_lmax = 8;
constexpr default_type default_msmt_normlambda = 1e-10;
constexpr default_type default_msmt_neglambda = 1e-10;

extern const App::OptionGroup MSMT_CSD_options;

class MSMT_CSD {
public:
  class Shared {
  public:
    Shared(const Header &dwi_header)
        : grad(DWI::get_DW_scheme(dwi_header)),
          HR_dirs(DWI::Directions::electrostatic_repulsion_300()),
          solution_min_norm_regularisation(default_msmt_normlambda),
          constraint_min_norm_regularisation(default_msmt_neglambda) {}

    void parse_cmdline_options() {
      using namespace App;
      auto opt = get_options("lmax");
      if (!opt.empty())
        lmax = parse_ints<uint32_t>(opt[0][0]);
      opt = get_options("directions");
      if (!opt.empty()) {
        const Eigen::MatrixXd directions = File::Matrix::load_matrix(opt[0][0]);
        DWI::Directions::validate(directions, opt[0][0], false);
        HR_dirs = Math::Sphere::as_spherical(directions);
      }
      opt = get_options("norm_lambda");
      if (!opt.empty())
        solution_min_norm_regularisation = opt[0][0];
      opt = get_options("neg_lambda");
      if (!opt.empty())
        constraint_min_norm_regularisation = opt[0][0];
      //TODO: parse GNL image
      opt = get_options("grad_nonlin");
      if (!opt.empty()){
        grad_dev = Image<float>::open(opt[0][0]);
        // check image dimensions are as expected
        if (grad_dev.ndim()!=4 || grad_dev.size(3)!=9)
          throw Exception("the gradient deviation image provided must be a 4D image with 9 volumes.");
        use_GNL = true;
      }
    }

    void set_responses(const std::vector<std::filesystem::path> &paths) {
      lmax_response.clear();
      for (const auto &p : paths) {
        Eigen::MatrixXd r;
        try {
          responses.push_back(File::Matrix::load_matrix(p));
        } catch (Exception &e) {
          try {
            se_responses.push_back(SEResponse(p)); // setsmodel param & the lmax for tissue
          } catch (Exception &e) {
            throw Exception(e, "File \"" + p.string() + "\" is not a valid response function file");
          }
        }
      }
      prepare_responses();
      response_files = paths;
    }

    // build the C matrix (option for voxel-dependent construction for GNL corrections)
    Eigen::MatrixXd build_C(size_t nparams,
                            uint32_t maxlmax,
                            const std::optional<Eigen::Vector3i>& vox = std::nullopt) {

      INFO("initialising multi-tissue CSD for " + str(num_tissues()) + " tissue types, with " + str(nparams) +
           " parameters");

      Eigen::MatrixXd C_local = Eigen::MatrixXd::Zero(grad.rows(), nparams);

      // gradient table
      Eigen::MatrixXd grad_local(grad.rows(), grad.cols());
      // in the case of GNL corrections:
      if (use_GNL && vox.has_value()) {
        GNL gnl(*this);
        Eigen::Matrix3d L_local;
        gnl.compute_L(L_local, *vox);
        gnl.correct_grad(grad_local, L_local);
      } else {
        // use a constant gradient table
        grad_local = grad;
      }

      std::vector<size_t> dwilist;
      for (size_t i = 0; i != static_cast<size_t>(grad.rows()); i++)
        dwilist.push_back(i);

      Eigen::MatrixXd directions = DWI::gen_direction_matrix(grad_local, dwilist);
      Eigen::MatrixXd SHT = Math::SH::init_transform(directions, maxlmax);
      for (ssize_t i = 0; i < SHT.rows(); i++)
        for (ssize_t j = 0; j < SHT.cols(); j++)
          if (std::isnan(SHT(i, j)))
            SHT(i, j) = 0.0;

      std::vector<size_t> shell_for_vol;

      if (responses.size()) {
        // assume shell structure:
        DWI::Shells shells(grad_local);
        shells.select_shells(false, false, false);

        for (size_t t = 0; t != num_tissues(); ++t) {
          if (static_cast<size_t>(responses[t].rows()) != shells.count())
            throw Exception("number of rows in response functions must match number of b-value shells; "
                            "number of shells is " +
                            str(shells.count()) + ", but file \"" + response_files[t].string() + "\" contains " +
                            str(responses[t].rows()) + " rows");
          // Pad response functions out to the requested lmax for this tissue
          responses[t].conservativeResizeLike(Eigen::MatrixXd::Zero(shells.count(), Math::ZSH::NforL(lmax[t])));
        }

        // TODO: is this just computing the Associated Legrendre polynomials...?
        // this is required for conversion from SH to RH:
        Eigen::MatrixXd delta(1, 2);
        delta << 0, 0;
        Eigen::VectorXd DSH_ = Math::SH::init_transform(delta, maxlmax).row(0);
        Eigen::VectorXd DSH(maxlmax / 2 + 1);
        size_t j = 0;
        for (ssize_t i = 0; i < DSH_.size(); i++)
          if (DSH_[i] != 0.0) {
            DSH[j] = DSH_[i];
            j++;
          }

        // convert responses from SH to RH:
        for (auto &r : responses)
          for (int c = 0; c < r.cols(); ++c)
            r.col(c) /= DSH[c];

        // reverse mapping from volume to shell index:
        shell_for_vol.resize(grad.rows());
        for (size_t shell_idx = 0; shell_idx < shells.count(); ++shell_idx) {
          const auto &vols = shells[shell_idx].get_volumes();
          for (size_t idx = 0; idx < vols.size(); idx++)
            shell_for_vol[vols[idx]] = shell_idx;
        }
      } else {
        // set up per-tissue SE responses:
        for (int t = 0; t != num_tissues(); ++t)
          se_responses[t].init(lmax[t]);
      }

      size_t pbegin = 0;
      for (size_t tissue_idx = 0; tissue_idx < num_tissues(); ++tissue_idx) {
        const size_t tissue_lmax = lmax[tissue_idx];
        const size_t tissue_n = Math::SH::NforL(tissue_lmax);
        const size_t tissue_nmzero = tissue_lmax / 2 + 1;
        Eigen::VectorXd fconv(tissue_n);
        Eigen::VectorXd workspace(tissue_n), se_R(tissue_n);

        for (size_t vol = 0; vol < grad.rows(); ++vol) {
          const size_t shell_idx = responses.size() ? shell_for_vol[vol] : 0;
          if (responses.empty())
            // computes se SH->RH coefficients for each (b,g)
            //TODO: if using GNL corrections - supply corrected bval for this voxel
            se_responses[tissue_idx].compute_SH_coeff(se_R, workspace, grad_local(vol,3));

          int li = 0;
          int mi = 0;
          for (int l = 0; l <= static_cast<int>(tissue_lmax); l += 2) {
            for (int m = -l; m <= l; m++) {
              // the rf matrix for shell-based structure
              fconv[mi] = responses.size() ? responses[tissue_idx](shell_idx, li) : se_R[li];
              mi++;
            }
            li++;
          }

          Eigen::VectorXd SHT_row(SHT.row(vol).head(tissue_n));
          SHT_row.array() *= fconv.array();
          C_local.row(vol).segment(pbegin, tissue_n) = SHT_row;
        }
        pbegin += tissue_n;
      }
      return C_local;
    }

    void init(const std::optional<Eigen::Vector3i>& vox = std::nullopt) {
      if (lmax.empty()) {
        lmax = lmax_response; // from prepare_response / SEResponse
        for (size_t t = 0; t != num_tissues(); ++t) {
          lmax[t] = std::min(default_msmt_lmax, lmax[t]);
        }
      } else {
        if (lmax.size() != num_tissues())
          throw Exception("Number of lmaxes specified (" + str(lmax.size()) +
              ") does not match number of tissues (" + str(num_tissues()) + ")");
        for (const auto i : lmax) {
          if (i % 2)
            throw Exception("Each value of lmax must be a non-negative even integer");
        }
      }
      // ensure GNL correction only implemented with SE responses
      if (use_GNL && responses.size())
        throw Exception("Gradient nonlinearity correction is only supported using the SE response function.");

      //////////////////////////////////////////////////
      // Set up the constrained least squares problem //
      //////////////////////////////////////////////////

      size_t nparams = 0;
      uint32_t maxlmax = 0;
      for (size_t i = 0; i < num_tissues(); i++) {
        nparams += Math::SH::NforL(lmax[i]);
        maxlmax = std::max(maxlmax, lmax[i]);
      }

      if (use_GNL && !vox.has_value())
        throw Exception("A voxel position is required to implement GNL-corrections.");
      Eigen::MatrixXd C;
      C = use_GNL ? build_C(nparams, maxlmax,vox) : build_C(nparams, maxlmax);

      // non-negativity constraint:
      std::vector<size_t> m(num_tissues());
      std::vector<size_t> n(num_tissues());
      size_t M = 0;
      size_t N = 0;

      Eigen::MatrixXd HR_SHT = Math::SH::init_transform(HR_dirs, maxlmax);

      for (size_t i = 0; i != num_tissues(); i++) {
        if (lmax[i] > 0)
          m[i] = HR_dirs.rows();
        else
          m[i] = 1;
        M += m[i];
        n[i] = Math::SH::NforL(lmax[i]);
        N += n[i];
      }

      Eigen::MatrixXd A(Eigen::MatrixXd::Zero(M, N));
      size_t b_m = 0;
      size_t b_n = 0;
      for (size_t i = 0; i != num_tissues(); i++) {
        A.block(b_m, b_n, m[i], n[i]) = HR_SHT.block(0, 0, m[i], n[i]);
        b_m += m[i];
        b_n += n[i];
      }

      //TODO: voxel-wise problem matrix for GNL corrections
      if (!use_GNL){
        problem = Math::ICLS::Problem<double>(
            C, A, Eigen::VectorXd(), 0, solution_min_norm_regularisation, constraint_min_norm_regularisation);

      INFO("Multi-shell, multi-tissue CSD initialised successfully");
      } else {
        // setup the problem for each voxel
      }
    }

    size_t num_tissues() const {
      if (responses.size())
        return responses.size();
      else if (se_responses.size())
        return se_responses.size();
      else
        throw Exception("no response defined");
    }

    const Eigen::MatrixXd grad;
    Eigen::MatrixXd HR_dirs;
    std::vector<uint32_t> lmax, lmax_response;
    std::vector<Eigen::MatrixXd> responses;
    std::vector<SEResponse> se_responses;
    std::vector<std::filesystem::path> response_files;
    Math::ICLS::Problem<double> problem;
    double solution_min_norm_regularisation, constraint_min_norm_regularisation;
    Image<float> grad_dev;
    bool use_GNL = false;

  private:
    void prepare_responses() {

      if (responses.size() && se_responses.size())
        throw Exception("cannot perform MSMT CSD using mixed response types");

      if (responses.size()) {
        for (size_t t = 0; t != num_tissues(); ++t) {
          Eigen::MatrixXd &r(responses[t]);
          size_t n = 0;
          for (Eigen::Index row = 0; row < r.rows(); row++) {
            for (Eigen::Index col = 0; col < r.cols(); col++) {
              if (r(row, col))
                n = std::max(n, static_cast<size_t>(col + 1));
            }
          }
          // Clip off any empty columns, i.e. degrees containing zero coefficients for all shells
          r.conservativeResize(r.rows(), n);
          // Store the lmax for each tissue based on their response functions;
          //   if the user doesn't manually specify lmax, these will determine the
          //   lmax of each tissue ODF output, with a further default lmax=8
          //   restriction at that stage
          lmax_response.push_back(Math::ZSH::LforN(r.cols()));
        }
      } else if (se_responses.size()) {
        for (const auto &r : se_responses)
          lmax_response.push_back(r.is_isotropic() ? 0 : default_msmt_lmax);
      } else
        throw Exception("no response defined");
    }
  };

  class GNL {
    public:
      GNL(const Shared &shared) : shared(shared) {}

      // compute the L tensor from grad_dev image at a voxel pos
      void compute_L(Eigen::Matrix3d &L, const Eigen::Vector3i& vox) const {
        // assuming L(x) stored as [Lxx, Lxy, Lxz, Lyx, Lyy, Lyz, Lzx, Lzy, Lzz]
        auto img = shared.grad_dev;

        img.index(0) = vox[0];
        img.index(1) = vox[1];
        img.index(2) = vox[2];

        for (int i = 0; i < 9; ++i) {
          img.index(3) = i;
          L(i / 3, i % 3) = img.value();
        }
      }

      // compute the corrected gradient information & bvalues at a voxel pos
      void correct_grad(Eigen::MatrixXd &grad_corr, const Eigen::Matrix3d &L) const {
        assert(grad_corr.size() == shared.grad.size());
        const Eigen::Matrix3d IL = Eigen::Matrix3d::Identity()+L;

        for (int N = 0; N < shared.grad.rows(); ++N) {
          const double b = shared.grad(N,3);
          const Eigen::Vector3d bv = shared.grad.row(N).head<3>();
          // account for left-handed system, i.e. flip the x-axis
          const Eigen::Vector3d v = flipMat() * bv ;
          const Eigen::Vector3d g = IL*v;
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
      const Shared &shared;

  };

  MSMT_CSD(const Shared &shared_data) : niter(0), shared(shared_data), solver(shared.problem) {}

  void operator()(const Eigen::VectorXd &data, Eigen::VectorXd &output) {
    if(!shared.use_GNL)
      niter = solver(output, data);
  }
  //TODO: voxel-wise solver solve

  size_t niter;
  const Shared &shared;

private:
  Math::ICLS::Solver<double> solver;
};

} // namespace MR::DWI::SDeconv

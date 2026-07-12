/*Copyright (c) 2008-2024 the MRtrix3 contributors.
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

#include <fstream>

#include "debug.h"
#include "file/matrix.h"
#include "math/ZSH.h"
#include "math/math.h"
#include "math/legendre.h"
#include "mrtrix.h"

namespace MR::DWI::SDeconv {

class SEResponse {

public:
  SEResponse() {}
  SEResponse(const std::string &filename) { load(filename); }
  SEResponse(SEResponse &&other) = default;
  //copy constructors for the gnl correction implementation (build C matrix)
  SEResponse(const SEResponse &other) = default;
  SEResponse &operator=(const SEResponse &other) = default;

  // for a given bval, compute the SH coefficients for RF
  void compute_SH_coeff(Eigen::VectorXd &sh_coeffs,Eigen::VectorXd &signals, const double bval) const {
    // se fit + ZSHT:
    if (is_isotropic()) {
      // param = [s0, D, alpha]
      double ex = (bval/1000.0) * se_coeffs[1];
      sh_coeffs(0) = se_coeffs[0] * exp(-std::pow(ex, se_coeffs[2]));
      //  dont scale by 1/sqrt(4*pi) - done by init_amp_transform
    } else {
      assert(signals.size() == sh_coeffs.size() && signals.size() == t_mat.rows());
      for (size_t i = 0; i < diffusivities.size(); i++)
        signals(i) = se_coeffs[0] * exp(-std::pow((bval/1000.0) * diffusivities[i], se_coeffs[3]));

      // convert amp (rf_signal) to zsh coefficients using iZSHT: coeffs = iZSHT * rf_signal
      sh_coeffs.noalias() = t_mat * signals;
    }

  }

  // check for isotropy (i.e. 3 model parameters)
  bool is_isotropic() const { return se_coeffs.size() == 3; }
 
  void init(int tissue_lmax) {
    if (tissue_lmax > 0 && is_isotropic())
      throw Exception("cannot use non-zero lmax for isotropic response");

    if (is_isotropic())
      return;

    size_t num_elev = Math::ZSH::NforL(tissue_lmax);
    diffusivities.resize(num_elev);
    std::vector<double> elevation(num_elev);

    Eigen::VectorXd alp (2*Math::ZSH::NforL(tissue_lmax));
    Math::Legendre::Plm_sph (alp, tissue_lmax, 0, 1.0);
    Eigen::VectorXd sh2rh (Math::ZSH::NforL(tissue_lmax));
    for (size_t l = 0; l <= tissue_lmax; l += 2)
      sh2rh[l/2] = 1.0 / alp[l];

    // evaluate at elevations evenly distributed between [0,pi/2]
    for (size_t i = 0; i < num_elev; i++){
      elevation[i] = (Math::pi / 2.0) * (double(i) / (num_elev - 1.0));
     // param = [s0, D_ax, D_rad, alpha]
      double cel = std::cos(elevation[i]);
      double sel = std::sin(elevation[i]);
      diffusivities[i] = se_coeffs[1] * (cel * cel) + se_coeffs[2] * (sel * sel);
    }
    // convert amp (rf_signal) to zsh coefficients using iZSHT, coeffs = iZSHT * rf_signal
    Eigen::MatrixXd transform = Math::ZSH::init_amp_transform<double>(elevation, tissue_lmax);
    // sh -> rh transformation
    t_mat = sh2rh.asDiagonal() * transform.inverse();

  }

  // load response file, set parameters and lmax:
  void load(const std::string &filename) {
    DEBUG("loading stretched exponential response from file \"" + filename + "\"");
    std::ifstream stream(filename, std::ios_base::binary);
    if (!stream)
      throw Exception("Unable to open stretched exponential coefficients file \"" + filename +
                      "\": " + strerror(errno));

    std::string magic;
    stream >> magic;
    if (magic != "SE")
      throw Exception("file \"" + filename + "\" is not in the expected format");

    double val;
    std::vector<double> coeffs; // the SE model parameters
    while (stream >> val)
      coeffs.push_back(val);

    if (coeffs.size() != 3 && coeffs.size() != 4)
      throw Exception("file \"" + filename +
                      "\" does not contain the expected number "
                      "of coefficients for a stretched exponential response");

    se_coeffs.resize(coeffs.size());
    for (size_t n = 0; n < coeffs.size(); ++n)
      se_coeffs[n] = coeffs[n];

    DEBUG("loaded stretched exponential response from file \"" + filename + "\" with coefficients " + str(coeffs));
  }

  size_t tissue_lmax;

private:
  Eigen::VectorXd se_coeffs, diffusivities;
  Eigen::MatrixXd t_mat; // transform matrix
};

} // namespace MR::DWI::SDeconv

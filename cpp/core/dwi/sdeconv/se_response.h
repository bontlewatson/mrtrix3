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
#include "mrtrix.h"

namespace MR::DWI::SDeconv {

class SEResponse {

public:
  SEResponse() {}
  SEResponse(const std::string &filename) { load(filename); }
  SEResponse(SEResponse &&other) = default;
  
  // for a given bval, compute the SH coefficients for RF 
  Eigen::VectorXd compute_SH_coeff(const double bval) const {
    // se fit + ZSHT: 
    if(is_isotropic()){
      // param = [s0, D, alpha]
      double ex = bval * se_coeffs[1];
      double rf_signal = se_coeffs[0] * exp(-std::pow(ex, se_coeffs[2]));
      //convert signal amp into ZSH coefficient (i.e. no directional dependence)
      // the coeff for lmax=0 is the signal amp?
      Eigen::VectorXd sh_coeffs(1);
      sh_coeffs(0) = rf_signal;
      return sh_coeffs;
    }
    else{
      // evaluate same signal at N number of directions (where N = N SH coefficinets)
      size_t num_elev = Math::ZSH::NforL(tissue_lmax);
      std::vector<double> elevations(num_elev);
      Eigen::VectorXd rf_signal(num_elev);

      // evaluate at elevations evenly distributed between [0,pi/2]
      for (size_t i = 0; i < num_elev; i++)
        elevations[i] = (Math::pi /2) * (double(i) / (num_elev - 1));

      // param = [s0, D_ax, D_rad, alpha]
      for (int i =0; i < num_elev; i++){
      double cel = std::cos(elevations[i]);
      double sel = std::sin(elevations[i]);
      double ex = bval * (se_coeffs[1]*(cel*cel) +se_coeffs[2]*(sel*sel));
      rf_signal(i) = se_coeffs[0]*exp(-std::pow(ex,se_coeffs[3]));
      }
      // convert amp (rf_signal) to zsh coefficients using iZSHT, coeffs = iZSHT * rf_signal
      Eigen::MatrixXd transform = Math::ZSH::init_amp_transform<double>(elevations, tissue_lmax);
      auto sh_coeffs = transform.colPivHouseholderQr().solve(rf_signal);
      return sh_coeffs;
    } 
  }

  // check for isotropy (i.e. 3 model parameters)
  bool is_isotropic() const { return se_coeffs.size() == 3; }
  // set an lmax according to no. of model parameters for tissue
  void set_lmax(){
    if (is_isotropic())
      tissue_lmax = 0;
    else
      tissue_lmax = 8; 
  }

  // load response file, set parameters and lmax: 
  void load (const std::string &filename) {
    DEBUG ("loading stretched exponential response from file \"" + filename + "\"");
    std::ifstream stream (filename, std::ios_base::binary);
    if (!stream)
      throw Exception("Unable to open stretched exponential coefficients file \""
          + filename + "\": " + strerror(errno));

    std::string magic;
    stream >> magic;
    if (magic != "SE")
      throw Exception ("file \"" + filename + "\" is not in the expected format");

    double val;
    std::vector<double> coeffs; // the SE model parameters
    while (stream >> val)
      coeffs.push_back (val);

    if (coeffs.size() != 3 && coeffs.size() != 4)
      throw Exception ("file \"" + filename + "\" does not contain the expected number "
          "of coefficients for a stretched exponential response");

    se_coeffs.resize (coeffs.size());
    for (size_t n = 0; n < coeffs.size(); ++n)
      se_coeffs[n] = coeffs[n];

    // once coefficients loaded, set tissue lmax 
    set_lmax();  

    DEBUG ("loaded stretched exponential response from file \"" + filename
        + "\" with coefficients " + str(coeffs));
  }
  
  size_t tissue_lmax;

private:
  Eigen::VectorXd se_coeffs;
  Eigen::VectorXd sh_coeffs;
};

}


/*
  Eigen::VectorXd coeffs(const double bval) {
    if (bval < original_bvals[0])
      throw Exception("bvalue out of bounds");

    if (bval >= original_bvals.back())
      return original_coeffs.row(original_coeffs.rows() - 1);

    size_t i = 0;
    while (bval > original_bvals[i + 1])
      i++;

    double ratio = (bval - original_bvals[i]) / (original_bvals[i + 1] - original_bvals[i]);
    return (1.0 - ratio) * original_coeffs.row(i) + ratio * original_coeffs.row(i + 1);
  }
*/
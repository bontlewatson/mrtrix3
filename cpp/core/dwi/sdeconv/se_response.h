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

  // TODO: init (lmax)
  // TODO: SH_coeffs (bval, coeffs&) const (!)

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
  bool is_isotropic() const { return se_coeffs.size() == 4; }

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
    std::vector<double> coeffs;
    while (stream >> val)
      coeffs.push_back (val);

    if (coeffs.size() != 3 && coeffs.size() != 4)
      throw Exception ("file \"" + filename + "\" does not contain the expected number "
          "of coefficients for a stretched exponential response");

    se_coeffs.resize (coeffs.size());
    for (size_t n = 0; n < coeffs.size(); ++n)
      se_coeffs[n] = coeffs[n];

    DEBUG ("loaded stretched exponential response from file \"" + filename
        + "\" with coefficients " + str(coeffs));
  }

private:
  Eigen::VectorXd se_coeffs;
};

}



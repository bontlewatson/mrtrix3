/* Copyright (c) 2008-2024 the MRtrix3 contributors.
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

#include "debug.h"
#include "file/matrix.h"
#include "math/ZSH.h"
#include "math/math.h"
#include "mrtrix.h"

namespace MR::DWI::SDeconv {

class Response {

public:
  Response() {}
  Response(const std::string &filename) { load(filename); }
  Response(Response &&other) = default;

  Eigen::VectorXd coeffs(const double bval) {
    if (bval < original_bvals[0])
      throw Exception("bvalue out of bounds");

    if (bval >= original_bvals.back())
      return original_coeffs.row(original_coeffs.rows() - 1);

    int i = 0;
    while (bval > original_bvals[i + 1])
      i++;

    float ratio = (bval - original_bvals[i]) / (original_bvals[i + 1] - original_bvals[i]);
    return (1.0 - ratio) * original_coeffs.row(i) + ratio * original_coeffs.row(i + 1);
  }

  int lmax() const { return Math::ZSH::LforN(original_coeffs.cols()); }

  const std::vector<double> &bvalues() const { return original_bvals; }

  // int size() const { return original_bvals.size(); }

  // LOAD MATRIX
  void load(const std::string &filename) {
    std::vector<std::string> comments;
    auto vec = File::Matrix::load_matrix_2D_vector(filename, &comments);

    // load bvalues
    for (const auto &comment : comments) {
      auto n = comment.rfind("Shells:");
      if (n != comment.npos)
        original_bvals = parse_floats(comment.substr(n + 7));
    }

    if (original_bvals.size() != vec.size())
      throw Exception("Number of b-values does not match the number of rows in response");

    if (!std::is_sorted(original_bvals.begin(), original_bvals.end()))
      throw Exception("B-values are not sorted in ascending order");

    original_coeffs.resize(vec.size(), vec[0].size());

    for (ssize_t r = 0; r < vec.size(); r++)
      for (ssize_t c = 0; c < vec[0].size(); c++)
        original_coeffs(r, c) = vec[r][c];
  }

private:
  Eigen::MatrixXd original_coeffs;
  std::vector<double> original_bvals;
};

} // namespace MR::DWI::SDeconv

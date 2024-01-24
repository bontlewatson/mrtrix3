/* Copyright (c) 2008-2023 the MRtrix3 contributors.
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

#ifndef __math_stats_types_h__
#define __math_stats_types_h__

#include "types.h"

#include <Eigen/Dense>

namespace MR
{
  namespace Math
  {
    namespace Stats
    {



      using index_type = uint32_t;
      using value_type = MR::default_type;
      using matrix_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
      using vector_type = Eigen::Array<value_type, Eigen::Dynamic, 1>;

      using mask_type = Eigen::Array<bool, Eigen::Dynamic, 1>;
      using index_array_type = Eigen::Array<index_type, Eigen::Dynamic, 1>;
      using shuffle_matrix_type = Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>;

      // Capability to internally store measurements at lower precision than
      //   that at which calculations are performed
      using measurements_value_type = float;
      using measurements_vector_type = Eigen::Matrix<measurements_value_type, Eigen::Dynamic, 1>;
      using measurements_matrix_type = Eigen::Matrix<measurements_value_type, Eigen::Dynamic, Eigen::Dynamic>;



    }
  }
}


#endif

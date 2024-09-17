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

#include "exception.h"
#include "file/mmap.h"
#include "file/utils.h"
#include <filesystem>

namespace MR::File {

inline void copy(const std::filesystem::path &source, const std::filesystem::path &destination) {
  {
    DEBUG("copying file \"" + source.string() + "\" to \"" + destination.string() + "\"...");
    MMap input(source);
    create(destination, input.size());
    MMap output(destination, true);
    ::memcpy(output.address(), input.address(), input.size());
  }
  check_app_exit_code();
}

} // namespace MR::File

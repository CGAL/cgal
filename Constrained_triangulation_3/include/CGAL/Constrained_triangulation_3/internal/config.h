// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CDT_3_CONFIG_H
#define CGAL_CDT_3_CONFIG_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/config.h>

#include <CGAL/Constrained_triangulation_3_types.h>

#if CGAL_CAN_USE_CXX20_FORMAT
#  define CGAL_CDT_3_CAN_USE_CXX20_FORMAT 1
#  include <format>
#endif

namespace CGAL {

constexpr bool cdt_3_msvc_2019_or_older() {
#if defined(_MSC_VER) && (_MSC_VER < 1930)
  return true;
#else
  return false;
#endif
}

#if CGAL_CDT_3_CAN_USE_CXX20_FORMAT

constexpr bool cdt_3_can_use_cxx20_format() {
  return true;
}

template <typename... Args>
decltype(auto) cdt_3_format(std::string_view fmt, const Args&... args) {
  return std::vformat(fmt, std::make_format_args(args...));
}

#else // not CGAL_CDT_3_CAN_USE_CXX20_FORMAT

template <typename... Args>
constexpr decltype(auto) cdt_3_format(Args&&...) {
  return std::string{};
}

constexpr bool cdt_3_can_use_cxx20_format() {
  return false;
}

#endif // not CGAL_CDT_3_CAN_USE_CXX20_FORMAT

} // namespace CGAL

#endif // CGAL_CDT_3_CONFIG_H

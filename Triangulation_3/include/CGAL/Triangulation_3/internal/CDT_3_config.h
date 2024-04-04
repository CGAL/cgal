// Copyright (c) 2019-2023  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CDT_3_CONFIG_H
#define CGAL_CDT_3_CONFIG_H

#include <CGAL/license/Triangulation_3.h>

#define CGAL_CDT_3_CAN_USE_CXX20_FORMAT 0
#if __has_include(<format>)
#  include <format>
#  if __cpp_lib_format >= 201907L || __cplusplus >= 202000L || _MSVC_LANG >= 202000L
#    undef CGAL_CDT_3_CAN_USE_CXX20_FORMAT
#    define CGAL_CDT_3_CAN_USE_CXX20_FORMAT 1
#  endif
#endif

namespace CGAL {
constexpr bool cdt_3_can_use_cxx20_format() {
#if CGAL_CDT_3_CAN_USE_CXX20_FORMAT
  return true;
#else
  return false;
#endif
}
} // namespace CGAL

#endif // CGAL_CDT_3_CONFIG_H

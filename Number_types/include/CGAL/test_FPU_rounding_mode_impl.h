// Copyright (c) 2008  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Laurent Rineau

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#ifndef CGAL_NDEBUG

#include <CGAL/basic.h>
#include <CGAL/assertions.h>

namespace CGAL {

struct Check_FPU_rounding_mode_is_restored {
  FPU_CW_t mode;

  Check_FPU_rounding_mode_is_restored()
    : mode( FPU_get_cw()) {}

  ~Check_FPU_rounding_mode_is_restored()
  {
    CGAL_assertion_msg( FPU_get_cw() == mode,
                        "The default FPU rounding mode has not been restored "
                        " before the exit of the program. "
                        "That may be a bug in some CGAL kernel code.");
  }
};

#ifdef CGAL_HEADER_ONLY

inline const Check_FPU_rounding_mode_is_restored& 
get_static_check_fpu_rounding_mode_is_restored()
{
  // A static object that emits a warning if the rounding mode at the
  // beginning and the end of the program are not the same.
  // Note that the get_static_check_fpu_rounding_mode_is_restored() function
  // must be called at least once so that this object is created.
  // It is done in the FPU_set_cw() function in FPU.h
  static const Check_FPU_rounding_mode_is_restored check_fpu_rounding_mode_is_restored;
  return check_fpu_rounding_mode_is_restored;
}

namespace {
  CGAL_UNUSED const Check_FPU_rounding_mode_is_restored &
    check_fpu_rounding_mode_is_restored 
    = get_static_check_fpu_rounding_mode_is_restored();
}

#else

// A global object that emits a warning if the rounding mode at the
// beginning and the end of the program are not the same.
static const Check_FPU_rounding_mode_is_restored check_fpu_rounding_mode_is_restored;

#endif // CGAL_HEADER_ONLY

} // end namespace CGAL

#endif // #ifnedef CGAL_NDEBUG

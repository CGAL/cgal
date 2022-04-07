// Copyright (c) 2008  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#ifndef CGAL_NDEBUG

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
                        "before the exit of the program. "
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

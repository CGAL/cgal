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
// 
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_NDEBUG

#include <CGAL/basic.h>
#include <CGAL/FPU.h>
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

// A global object that emits a warning if the rounding mode at the
// beginning and the end of the program are not the same.
static const Check_FPU_rounding_mode_is_restored check_fpu_rounding_mode_is_restored;

} // end namespace CGAL

#endif // #ifnedef CGAL_NDEBUG

// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// 
//
// Author(s)     : Kaspar Fischer


#ifndef CGAL_MINIBALL_CONFIGURE
#define CGAL_MINIBALL_CONFIGURE

// Remark: In case you want to fine-tune the code, feel free to change
// the options below.

// Option: Namespace name
//
// Background: By default, all data-structures and routines of the
// package are placed in a namespace called CGAL.  You can
// change this name by altering the following #define.
//
// Default value: CGAL
#define CGAL_MINIBALL_NAMESPACE CGAL

// Option: Assertions
//
// Background: The package contains lots of assertions (i.e. internal
// consistency checks).  For instance, if assertions are enabled then
// the package will complain when you add balls with negative radii.
// If you disable assertions, such tests will not be made (with the
// advantage that the code is slightly faster).  Do *not* disable
// assertions during development.
//
// Default setting: defined
#ifndef CGAL_NO_ASSERTIONS
#define CGAL_MINIBALL_DEBUG
#endif
#ifdef NDEBUG
#undef  CGAL_MINIBALL_DEBUG
#endif

// (You should not have to alter anything below here.)

// If CGAL is not being used, we need to define certain things:
#ifndef CGAL_VERSION
  namespace CGAL_MINIBALL_NAMESPACE {
    struct Tag_true {};
    struct Tag_false {};
  }
  #define CGAL_MINIBALL_NTS 
#else
  #include <CGAL/basic.h>
  #define CGAL_MINIBALL_NTS CGAL_NTS
#endif

// Define some assertion macros used in the code.
#ifdef CGAL_MINIBALL_DEBUG
  #include <CGAL/assertions.h>
  #define CGAL_MINIBALL_ASSERT(expr) CGAL_assertion(expr)
  #define CGAL_MINIBALL_DO_DEBUG(expr) expr
#else
  #define CGAL_MINIBALL_ASSERT(expr) ;
  #define CGAL_MINIBALL_DO_DEBUG(expr) ;
#endif

// Currently, we include all code in the header files because most
// compilers don't support exporting templates anyway:
#define CGAL_MINIBALL_NO_TEMPLATE_EXPORT

#include <bitset>
#include <sstream>

#endif // CGAL_MINIBALL_CONFIGURE

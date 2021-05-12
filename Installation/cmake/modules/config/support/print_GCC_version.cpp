// Copyright (c) 2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Fernando Cacciola

// Print out the gcc version

#include <iostream>


//
// Just in case this is called with a non-gcc compiler such as pgCC
//

#ifdef __clang_version__
#  define _CGAL_GCC_VERSION "Not GNU/CC (CLANG)"
#else
#  ifndef __GNUC__
#    define _CGAL_GCC_VERSION "Not GNU/CC"
#  endif
#endif
#ifndef _CGAL_GCC_VERSION
#  ifdef __VERSION__
#    define _CGAL_GCC_VERSION __VERSION__
#  else
#    define _CGAL_GCC_VERSION "Unknown version (_CGAL_GCC_VERSION is not defined)"
#  endif
#endif

int main()
{
  std::cout << "version=" << _CGAL_GCC_VERSION << std::endl;
  return 0;
}

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

#ifndef __INTEL_COMPILER
#  define CGAL_INTEL_VERSION "Not ICC"
#else
#  define CGAL_INTEL_VERSION __INTEL_COMPILER
#endif

int main()
{
  std::cout << "version=" << CGAL_INTEL_VERSION << std::endl;
  return 0;
}

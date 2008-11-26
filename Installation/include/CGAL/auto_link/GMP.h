// Copyright (c) 2007-2008 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Fernando Cacciola (fernando.cacciola@geometryfactry.com)

#ifndef CGAL_AUTO_LINK_GMP_H
#define CGAL_AUTO_LINK_GMP_H

#ifndef CGAL_NO_AUTOLINK_GMP

#define CGAL_LIB_NAME gmp

// One hardcodes a toolset for gmp. We do not need to distinguish between
// vc8 or vc9 (or even vc10), because gmp is a C library (not C++).
#define CGAL_LIB_TOOLSET vc80

#include <CGAL/auto_link/auto_link.h>

#endif // CGAL_NO_AUTOLINK_GMP

#endif // CGAL_AUTO_LINK_GMP_H


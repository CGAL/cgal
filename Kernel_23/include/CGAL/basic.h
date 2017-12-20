// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Lutz Kettner
//                 Stefan Schirra
 
#ifndef CGAL_BASIC_H
#define CGAL_BASIC_H

#include <CGAL/config.h>

#include <iostream>
#include <cstdlib>

// This cannot be disabled for now until we have a clear idea which
// compilers implement N3276.

// #if !defined(CGAL_CFG_NO_CPP0X_DECLTYPE)
//   #define BOOST_RESULT_OF_USE_DECLTYPE
// #endif
#include <CGAL/result_of.h>

#include <CGAL/assertions.h>
#include <CGAL/tags.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/IO/io.h>
#include <CGAL/kernel_basic.h>

#endif // CGAL_BASIC_H

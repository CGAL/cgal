// Copyright (c) 2011
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
// Author     : Laurent Rineau

#ifndef CGAL_VERSION_MACROS_H
#define CGAL_VERSION_MACROS_H

#include <CGAL/version.h>

#ifndef CGAL_STR
#define CGAL_STR(X) CGAL_STR_STR(X)
#define CGAL_STR_STR(X) #X
#endif

#ifndef CGAL_str
#define CGAL_xstr(s) #s
#define CGAL_str(s) CGAL_xstr(s)
#endif

#define CGAL_VERSION_STR CGAL_STR(CGAL_VERSION)

// The following macro definitions:
//   - do not use extra parenthesis,
//   - and do not use whitespace 
// on purpose, so that the Windows Resource Compiler can understand 
// the file generated from src/CGAL_libs_verinfo.rc.in
#define CGAL_VERSION_MAJOR (CGAL_VERSION_NR/10000000%100)
#define CGAL_VERSION_MINOR (CGAL_VERSION_NR/100000%100)
#define CGAL_VERSION_PATCH (CGAL_VERSION_NR/10000%10)
#define CGAL_VERSION_BUILD (CGAL_VERSION_NR%10000)

#endif // CGAL_VERSION_MACROS_H

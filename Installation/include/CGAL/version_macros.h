// Copyright (c) 2011
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

// Copyright (c) 2007 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Fernando Cacciola (fernando.cacciola@geometryfactry.com)

#ifndef CGAL_AUTO_LINK_CGAL_H
#define CGAL_AUTO_LINK_CGAL_H

#include <CGAL/config.h>

#ifndef CGAL_NO_AUTOLINK_CGAL
#ifndef CGAL_EXPORTS
// If CGAL_EXPORTS is defined it means that we are building the CGAL
// library as a DLL.

#define CGAL_LIB_NAME CGAL
#include <CGAL/auto_link/auto_link.h>

#endif // CGAL_EXPORTS
#endif // CGAL_NO_AUTOLINK_CGAL

#endif // CGAL_AUTO_LINK_CGAL_H

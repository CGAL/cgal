// Copyright (c) 2016  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_MSVC_COMPILER_CONFIG_H
#define CGAL_MSVC_COMPILER_CONFIG_H

// For all known version of MSVC. Actually we do not really have a
// test for that bug.
#define CGAL_CFG_MATCHING_BUG_6 1

// Fixed since MSVC 2015
#define CGAL_CFG_MATCHING_BUG_7 1

// for all known version of MSVC
#define CGAL_CFG_MATCHING_BUG_8 1

// Should be only for MSVC 2012 and 2013 in Release
#define CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG 1

#endif // CGAL_MSVC_COMPILER_CONFIG_H

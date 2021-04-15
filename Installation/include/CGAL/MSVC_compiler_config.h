// Copyright (c) 2016  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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

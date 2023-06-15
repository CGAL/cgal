// Copyright (c) 2023 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : -

#ifndef CGAL_VERSION_ENFORCER_H
#define CGAL_VERSION_ENFORCER_H

#include <CGAL/version_macros.h>

// All files including this header are meant to work with a given version of CGAL
// When using forked headers, set the 4 following macros to the version of CGAL
// you want to use.
#define CGAL_AUTHORIZED_VERSION CGAL_VERSION_STR
#define CGAL_AUTHORIZED_VERSION_MAJOR CGAL_VERSION_MAJOR
#define CGAL_AUTHORIZED_VERSION_MINOR CGAL_VERSION_MINOR
#define CGAL_AUTHORIZED_VERSION_PATCH CGAL_VERSION_PATCH

// Check that the version of CGAL used is the one expected
#if (CGAL_VERSION_MAJOR != CGAL_AUTHORIZED_VERSION_MAJOR)
#pragma message "You are using CGAL version: " CGAL_STR(CGAL_VERSION) "."
#error This header is meant to be with used with CGAL 5.5 only.
#endif

#if (CGAL_VERSION_MINOR != CGAL_AUTHORIZED_VERSION_MINOR)
#pragma message "You are using CGAL version: " CGAL_STR(CGAL_VERSION) "."
#error This header is meant to be with used with CGAL 5.5 only.
#endif

#if (CGAL_VERSION_PATCH != CGAL_AUTHORIZED_VERSION_PATCH)
#pragma message "You are using CGAL version: " CGAL_STR(CGAL_VERSION) "."
#error This header is meant to be with used with CGAL 5.5 only.
#endif

#endif // CGAL_VERSION_ENFORCER_H

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

#ifndef CGAL_VERSION_CHECKER_H
#define CGAL_VERSION_CHECKER_H

#include <CGAL/version_macros.h>

// All files including this header are meant to work with a given version of CGAL
// When using forked headers, set the following macro to the version of CGAL
// you want to use.

//// Set the 3 following macros to the version of CGAL you want to use
//#define CGAL_COMPATIBLE_VERSION_MAJOR 6
//#define CGAL_COMPATIBLE_VERSION_MINOR 0
//#define CGAL_COMPATIBLE_VERSION_PATCH 0

// Set the following macros to 1 to get a warning/an error
// when using a possibly incompatible version of CGAL
#define CGAL_VERSION_CHECKER_ERROR 0
#define CGAL_VERSION_CHECKER_WARNING 0

#define CGAL_COMPATIBLE_VERSION_STR CGAL_STR(CGAL_COMPATIBLE_VERSION_MAJOR) "." \
                                    CGAL_STR(CGAL_COMPATIBLE_VERSION_MINOR) "." \
                                    CGAL_STR(CGAL_COMPATIBLE_VERSION_PATCH)


// Check that the version of CGAL used is the one expected
#if  CGAL_COMPATIBLE_VERSION_MAJOR != CGAL_VERSION_MAJOR \
  || CGAL_COMPATIBLE_VERSION_MINOR != CGAL_VERSION_MINOR \
  || CGAL_COMPATIBLE_VERSION_PATCH != CGAL_VERSION_PATCH

    #if CGAL_VERSION_CHECKER_WARNING || CGAL_VERSION_CHECKER_ERROR
      #pragma message("These headers are meant to be used with CGAL " CGAL_COMPATIBLE_VERSION_STR " only."\
                       " You are using CGAL version: " CGAL_STR(CGAL_VERSION) ".")

      #ifdef CGAL_VERSION_CHECKER_ERROR
        #error "Incompatible version of CGAL"
      #endif

    #endif

#endif

#endif // CGAL_VERSION_CHECKER_H

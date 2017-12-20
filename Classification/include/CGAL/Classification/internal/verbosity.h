// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CLASSIFICATION_INTERNAL_VERBOSITY_H
#define CLASSIFICATION_INTERNAL_VERBOSITY_H

#include <CGAL/license/Classification.h>

// General verbosity

#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_SILENT false
#else
#define CGAL_CLASSIFICATION_SILENT true
#endif

#define CGAL_CLASSIFICATION_CERR \
  if(CGAL_CLASSIFICATION_SILENT) {} else std::cerr

// Verbosity for training part

#if defined(CGAL_CLASSTRAINING_VERBOSE)
#define CGAL_CLASSTRAINING_SILENT false
#else
#define CGAL_CLASSTRAINING_SILENT true
#endif

#define CGAL_CLASSTRAINING_CERR \
  if(CGAL_CLASSTRAINING_SILENT) {} else std::cerr

#endif

// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

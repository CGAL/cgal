// Copyright (c) 1999,2007,2012
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
// Author(s)     : Stefan Schirra, Michael Hemmer

#ifndef CGAL_NUMBER_TYPE_CONFIG_H
#define CGAL_NUMBER_TYPE_CONFIG_H

#include <CGAL/config.h>

#define CGAL_PI 3.141592653589793238462643383279502884
#define CGAL_SQRT2 1.414213562373095048801688724209698078
#define CGAL_SQRT3 1.732050807568877293527446341505872366
#define CGAL_SQRT5 2.236067977499789696409173668731276235


#ifdef CGAL_USE_NTS_NAMESPACE

#define CGAL_NTS_BEGIN_NAMESPACE namespace NTS {
#define CGAL_NTS_END_NAMESPACE }
#define CGAL_NTS ::CGAL::NTS::

#else

#define CGAL_NTS_BEGIN_NAMESPACE
#define CGAL_NTS_END_NAMESPACE
#define CGAL_NTS ::CGAL::

#endif

#endif // CGAL_NUMBER_TYPE_CONFIG_H

// Copyright (c) 1997-2001
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_DEBUG_H
#define CGAL_OPTIMISATION_DEBUG_H

// macro definitions
// =================

// debug
// -----
#if (    defined( CGAL_OPTIMISATION_NO_DEBUG) \
      || defined( CGAL_NO_DEGUG) || defined( NDEBUG))
#  define  CGAL_optimisation_debug  if ( 0)
#else
#  define  CGAL_optimisation_debug  if ( 1)
#endif // optimisation debug

#endif // CGAL_OPTIMISATION_DEBUG_H

// ===== EOF ==================================================================

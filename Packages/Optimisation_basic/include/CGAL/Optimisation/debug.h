// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Optimisation/debug.h
// package       : Optimisation_basic 3.8.3 (20 Mar 2001)
// chapter       : Geometric Optimisation
//
// source        : web/Optimisation_basic.aw
// revision      : 5.9
// revision_date : 2001/03/20 16:54:14
//
// author(s)     : Sven Schönherr
// maintainer    : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: debug macro for optimisation algorithms
// ============================================================================

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

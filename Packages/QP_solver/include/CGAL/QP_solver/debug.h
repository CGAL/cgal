// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// implementation: debug macro for the QP engine
// ======================================================================

#ifndef CGAL_QPE_DEBUG_H
#define CGAL_QPE_DEBUG_H

// macro definitions
// =================

// debug
// -----
#if (    defined( CGAL_QPE_NO_DEBUG) \
      || defined( CGAL_NO_DEGUG) || defined( NDEBUG))
#  define  CGAL_qpe_debug  if ( 0)
#else
#  define  CGAL_qpe_debug  if ( 1)
#endif // qpe debug

#endif // CGAL_QPE_DEBUG_H

// ===== EOF ==================================================================

// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/LEDA_basic.h
// package       : LEDA (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0
// revision_date : 19 March 2002
// author(s)     : Matthias Baesken
//                
//
// coordinator   :  Matthias Baesken <baesken@informatik.uni-trier.de>
// ======================================================================

#ifndef CGAL_LEDA_BASIC_H
#define CGAL_LEDA_BASIC_H

#ifdef CGAL_USE_LEDA
// The following is needed for LEDA 4.4 due to min/max problems...
#  define LEDA_NO_MIN_MAX_TEMPL

#include <LEDA/basic.h>

#ifdef LEDA_NAMESPACE
#  define CGAL_LEDA_SCOPE  leda
#else
#  define CGAL_LEDA_SCOPE 
#endif


#endif


#endif

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
// release_date  : 2001, May 23
//
// file          : include/CGAL/LEDA/global.h
// package       : cgal_window (0.9.7)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_GLOBAL_H
#define CGAL_WINDOW_GLOBAL_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

//------------------------------------------------------------------------------
// global types, constants, and macros
//------------------------------------------------------------------------------

typedef void* GenPtr;    // generic pointer type

#ifndef nil
#define nil 0
#endif


#define	LEDA_PI   3.14159265358979323846
#define	LEDA_PI_2 1.57079632679489661923

//------------------------------------------------------------------------------
// values
//------------------------------------------------------------------------------

#if defined(__USE_VALUES_H__)

#include <values.h>

#else

#include <cfloat>
#include <climits>


#if !defined(MAXINT)
#define MAXINT INT_MAX
#endif

#if !defined(MAXDOUBLE)
#define MAXDOUBLE DBL_MAX
#endif

#endif


//------------------------------------------------------------------------------
// LEDA class 
//------------------------------------------------------------------------------

namespace CGAL {

struct __exportC LEDA {

static const char* copyright_string;
static const char* copyright_window_string;
};

}

#endif

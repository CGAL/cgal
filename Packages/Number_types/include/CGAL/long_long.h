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
// release_date  : 2000, December 13
// 
// source        : webNT/Long_long.fw
// file          : include/CGAL/long_long.h
// package       : Number_types (4.2)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 4.2
// revision_date : 13 Dec 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

// !!!  ANSI C++ does not support `long long'  !!!

#ifndef CGAL_LONG_LONG_H
#define CGAL_LONG_LONG_H

#include <CGAL/number_type_tags.h>

CGAL_BEGIN_NAMESPACE


inline
double
to_double(long long int i)
{ return (double)i; }

inline
Number_tag
number_type_tag(long long int)
{ return Number_tag(); }

inline
bool
is_finite(long long int)
{ return true; }

inline
bool
is_valid(long long int)
{ return true; }

#if (defined(__sparc__) || defined(__sparc) || defined(sparc)) || \
    (defined(__sgi__)   || defined(__sgi)   || defined(sgi)) || \
    (defined(__i386__)  || defined(__i386)  || defined(i386)) || \
    (defined(__powerpc__) || defined(__powerpc) || defined(powerpc))
typedef  long long int           Integer64;
typedef  unsigned long long int  UInteger64;
#endif // __sparc__ || ...

CGAL_END_NAMESPACE

#endif // CGAL_LONG_LONG_H

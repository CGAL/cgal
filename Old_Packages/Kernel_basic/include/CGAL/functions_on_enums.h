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
// release_date  : 2000, December 10
// 
// source        : functions_on_enums.fw
// file          : functions_on_enums.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.17
// revision_date : 10 Dec 2000 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_FUNCTIONS_ON_ENUMS_H
#define CGAL_FUNCTIONS_ON_ENUMS_H

#ifndef CGAL_CONFIG_H
#include <CGAL/config.h>
#endif  // CGAL_CONFIG_H

#ifndef CGAL_ENUM_H
#include <CGAL/enum.h>
#endif // CGAL_ENUM_H

CGAL_BEGIN_NAMESPACE

template <class T>
inline
T
opposite(const T& t)
{ return -t; }

CGAL_TEMPLATE_NULL
inline
Sign
opposite(const Sign& o)
{ return (Sign)(-(int)o); }

CGAL_TEMPLATE_NULL
inline
Oriented_side
opposite(const Oriented_side& os)
{ return (Oriented_side)(-(int)os); }

CGAL_TEMPLATE_NULL
inline
Bounded_side
opposite(const Bounded_side &bs)
{ return (Bounded_side)(-(int)bs); }

CGAL_END_NAMESPACE


#endif // CGAL_FUNCTIONS_ON_ENUMS_H

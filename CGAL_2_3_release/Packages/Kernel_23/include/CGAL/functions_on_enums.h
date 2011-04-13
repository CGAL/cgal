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
// file          : functions_on_enums.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_FUNCTIONS_ON_ENUMS_H
#define CGAL_FUNCTIONS_ON_ENUMS_H

#include <CGAL/config.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template <class T>
inline
T
opposite(const T& t)
{ return -t; }

inline
Sign
opposite(Sign o)
{ return static_cast<Sign>( - static_cast<int>(o)); }

inline
Oriented_side
opposite(Oriented_side os)
{ return static_cast<Oriented_side>( - static_cast<int>(os)); }

inline
Bounded_side
opposite(Bounded_side bs)
{ return static_cast<Bounded_side>( - static_cast<int>(bs)); }

inline
Angle
opposite(Angle a)
{ return static_cast<Angle>( - static_cast<int>(a)); }

CGAL_END_NAMESPACE

#endif // CGAL_FUNCTIONS_ON_ENUMS_H

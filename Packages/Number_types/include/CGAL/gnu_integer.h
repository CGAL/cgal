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
// file          : gnu_integer.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_GNU_INTEGER_H
#define CGAL_GNU_INTEGER_H

#ifndef CGAL_PROTECT_INTEGER_H
#include <Integer.h>
#define CGAL_PROTECT_INTEGER_H
#endif // CGAL_PROTECT_INTEGER_H

CGAL_BEGIN_NAMESPACE

inline
bool
is_finite(const Integer &)
{ return true; }

inline
bool
is_valid(const Integer &)
{ return true; }

inline
double
to_double(const Integer & i)
{ return i.as_double(); }

inline
Number_tag
number_type_tag(const Integer& )
{ return Number_tag(); }

inline
io_Operator
io_tag(const Integer& )
{ return io_Operator(); }

CGAL_END_NAMESPACE

#endif // CGAL_GNU_INTEGER_H

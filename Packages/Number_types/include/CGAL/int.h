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
// source        : Int.fw
// file          : int.h
// package       : Number_types (4.2)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 4.2
// revision_date : 13 Dec 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_INT_H
#define CGAL_INT_H

#ifndef CGAL_NUMBER_TYPE_TAGS_H
#include <CGAL/number_type_tags.h>
#endif // CGAL_NUMBER_TYPE_TAGS_H

// int

CGAL_BEGIN_NAMESPACE


inline
double
to_double(int i)
{ return (double)i; }

inline
Interval_base
to_interval(int i)
{ return CGAL::to_interval(double(i)); }

inline
Number_tag
number_type_tag(int)
{ return Number_tag(); }

inline
bool
is_finite(int)
{ return true; }

inline
bool
is_valid(int)
{ return true; }

inline
io_Operator
io_tag(int)
{ return io_Operator(); }

// long

inline
double
to_double(long int i)
{ return (double)i; }

inline
Interval_base
to_interval(long int i)
{
  // actually we would like to compare number of mantissa bits,
  // this seems to be a sufficient approximation
  CGAL_kernel_assertion( sizeof( double) > sizeof( long));
  // need something else for 64 bits longs.
  return CGAL::to_interval(double(i));
}

inline
Number_tag
number_type_tag(long int)
{ return Number_tag(); }

inline
bool
is_finite(long int)
{ return true; }

inline
bool
is_valid(long int)
{ return true; }

inline
io_Operator
io_tag(long int)
{ return io_Operator(); }

// short

inline
double
to_double(short int i)
{ return (double)i; }

inline
Interval_base
to_interval(short int i)
{ return CGAL::to_interval(double(i)); }

inline
Number_tag
number_type_tag(short int)
{ return Number_tag(); }

inline
bool
is_finite(short int)
{ return true; }

inline
bool
is_valid(short int)
{ return true; }

inline
io_Operator
io_tag(short int)
{ return io_Operator(); }

// long long

#ifdef LONG_LONG

inline
double
to_double(long long i)
{ return (double)i; }

inline
Interval_base
to_interval(long int i)
{
  // actually we would like to compare number of mantissa bits,
  // this seems to be a sufficient approximation
  CGAL_kernel_assertion( sizeof( double) > sizeof( long long));
  // need something else for 64 bits longs.
  return CGAL::to_interval(double(i));
}

inline
Number_tag
number_type_tag(long long)
{ return Number_tag(); }

inline
bool
is_finite(long long)
{ return true; }

inline
bool
is_valid(long long)
{ return true; }
#endif // LONG_LONG

// io_tags for unsigned types
inline
io_Operator
io_tag(unsigned char)
{ return io_Operator(); }

inline
io_Operator
io_tag(unsigned short)
{ return io_Operator(); }

inline
io_Operator
io_tag(unsigned int)
{ return io_Operator(); }

inline
io_Operator
io_tag(unsigned long)
{ return io_Operator(); }

CGAL_END_NAMESPACE


#endif // CGAL_INT_H

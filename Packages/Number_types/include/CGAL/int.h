// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_INT_H
#define CGAL_INT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

// int

template <> struct Number_type_traits<int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

inline
int
div(int i1, int i2)
{ return i1 / i2; }

inline
double
to_double(int i)
{ return static_cast<double>(i); }

inline 
std::pair<double,double>
to_interval(int i)
{
  return std::pair<double,double>(i, i);
}

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

template <> struct Number_type_traits<long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

inline
long int
div(long int i1, long int i2)
{ return i1 / i2; }

inline
double
to_double(long int i)
{ return static_cast<double>(i); }

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

template <> struct Number_type_traits<short int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

inline
short int
div(short int i1, short int i2)
{ return i1 / i2; }

inline
double
to_double(short int i)
{ return static_cast<double>(i); }

inline 
std::pair<double,double>
to_interval(short int i)
{
  return std::pair<double,double>(i, i);
}

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


// Note : "long long" support is in <CGAL/long_long.h>

CGAL_END_NAMESPACE

#endif // CGAL_INT_H

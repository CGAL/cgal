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
#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

// int

#ifdef CGAL_NEW_NT_TRAITS

template <>
struct Number_type_traits<int>
  : public CGALi::Default_euclidean_ring_number_type_traits<int>
{
  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  typedef Tag_false  Has_simplify;

  static inline double to_double(int i) {
    return static_cast<double>(i);
  }

  static inline std::pair<double,double> to_interval(int i) {
    // this may not be correct since an integer may not be
    // representable by a double
    return std::pair<double,double>(i, i);
  }

  static inline bool is_finite(int) { return true; }
  static inline bool is_valid(int)  { return true; }
};

#else

template <> struct Number_type_traits<int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
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

// unsigned int

template <> struct Number_type_traits<unsigned int> {
  typedef Tag_true  Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

inline
unsigned int
div(unsigned int i1, unsigned int i2)
{ return i1 / i2; }

inline
double
to_double(unsigned int i)
{ return static_cast<double>(i); }

inline
bool
is_finite(unsigned int)
{ return true; }

inline
bool
is_valid(unsigned int)
{ return true; }

inline
io_Operator
io_tag(unsigned int)
{ return io_Operator(); }

inline
unsigned int
is_negative(unsigned int)
{ return false; }
  
inline
Sign
sign(unsigned int i) 
{ return is_positive(i) ? POSITIVE : ZERO; }
  
inline
unsigned int
is_positive(unsigned int i)
{ return i != 0; }

inline
unsigned int
abs(unsigned int i)
{ return i; }

#endif

// long

#ifdef CGAL_NEW_NT_TRAITS

template <>
struct Number_type_traits<long int>
  : public CGALi::Default_euclidean_ring_number_type_traits<long int>
{
  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  static inline double to_double(long int i) {
    return static_cast<double>(i);
  }

  static inline bool is_finite(long int i) { return true; }
  static inline bool is_valid(long int i)  { return true; }

  static inline std::pair<double,double>
  to_interval(const long int& l) {
#ifndef __BORLANDC__ // The stupid Borland compiler generates warnings...
    if (sizeof(double) > sizeof(long)) {
      // On 64bit platforms, a long doesn't fit exactly in a double.
      // Well, a perfect fix would be to use std::numeric_limits<>, but...
      Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
      Interval_nt<false> approx ((double) l);
      FPU_set_cw(CGAL_FE_UPWARD);
      approx += Interval_nt<false>::smallest();
      return approx.pair();
    }
    else
#endif
      return std::pair<double,double>(l,l);
  }
};

#else // CGAL_NEW_NT_TRAITS

template <> struct Number_type_traits<long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
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

// long

template <> struct Number_type_traits<unsigned long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

inline
unsigned long int
div(unsigned long int i1, unsigned long int i2)
{ return i1 / i2; }

inline
double
to_double(unsigned long int i)
{ return static_cast<double>(i); }

inline
bool
is_finite(unsigned long int)
{ return true; }

inline
bool
is_valid(unsigned long int)
{ return true; }

inline
io_Operator
io_tag(unsigned long int)
{ return io_Operator(); }

inline
unsigned long int
is_negative(unsigned long int) 
{ return false; }
  
inline
Sign
sign(unsigned long int i) 
{ return is_positive(i) ? POSITIVE : ZERO; }
  
inline
unsigned long int
is_positive(unsigned long int i)
{ return i != 0; }

inline
unsigned long int
abs(unsigned long int i)
{ return i; }

#endif // CGAL_NEW_NT_TRAITS

// short

#ifdef CGAL_NEW_NT_TRAITS

template<>
struct Number_type_traits<short int>
  : public CGALi::Default_euclidean_ring_number_type_traits<short int>
{
  typedef Tag_false Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;

  typedef Tag_false Has_simplify;

  static inline double to_double(short int i) {
    return static_cast<double>(i);
  }

  static inline std::pair<double,double> to_interval(short int i) {
    return std::pair<double,double>(i, i);
  }

  static inline bool is_finite(short int) { return true; }
  static inline bool is_valid(short int)  { return true; }
};

#else // CGAL_NEW_NT_TRAITS

template <> struct Number_type_traits<short int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
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

// unsigned short

template <> struct Number_type_traits<unsigned short int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

inline
unsigned short int
div(unsigned short int i1, unsigned short int i2)
{ return i1 / i2; }

inline
double
to_double(unsigned short int i)
{ return static_cast<double>(i); }

inline
bool
is_finite(unsigned short int)
{ return true; }

inline
bool
is_valid(unsigned short int)
{ return true; }

inline
io_Operator
io_tag(unsigned short int)
{ return io_Operator(); }

inline
unsigned short int
is_negative(unsigned short int) 
{ return false; }
  
inline
Sign
sign(unsigned short int i) 
{ return is_positive(i) ? POSITIVE : ZERO; }
  
inline
unsigned short int
is_positive(unsigned short int i)
{ return i != 0; }

inline
unsigned short int abs(unsigned short int i)
{ return i; }

#endif // CGAL_NEW_NT_TRAITS

// Note : "long long" support is in <CGAL/long_long.h>

CGAL_END_NAMESPACE

#endif // CGAL_INT_H

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
// file          : int.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_INT_H
#define CGAL_INT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

// int

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

namespace NTS {
  inline unsigned int is_negative(unsigned int i) { return false; }
  
  inline unsigned int sign(unsigned int i) 
    { return is_positive(i) ? POSITIVE : ZERO; }
  
  inline unsigned int is_positive(unsigned int i)
    { return i != 0; }

  inline unsigned int abs(unsigned int i) { return i; }
} // namespace NTS

// long

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

namespace NTS {
  inline unsigned long int is_negative(unsigned long int i) 
    { return false; }
  
  inline unsigned long int sign(unsigned long int i) 
    { return is_positive(i) ? POSITIVE : ZERO; }
  
  inline unsigned long int is_positive(unsigned long int i)
    { return i != 0; }

  inline unsigned long int abs(unsigned long int i) { return i; }
} // namespace NTS

// short

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

namespace NTS {
  inline unsigned short int is_negative(unsigned short int i) 
    { return false; }
  
  inline unsigned short int sign(unsigned short int i) 
    { return is_positive(i) ? POSITIVE : ZERO; }
  
  inline unsigned short int is_positive(unsigned short int i)
    { return i != 0; }

  inline unsigned short int abs(unsigned short int i) { return i; }
} // namespace NTS

#ifdef LONG_LONG

// long long

template <> struct Number_type_traits<long long> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

inline
long long
div(long long i1, long long i2)
{ return i1 / i2; }

inline
double
to_double(long long i)
{ return static_cast<double>(i); }

inline
bool
is_finite(long long)
{ return true; }

inline
bool
is_valid(long long)
{ return true; }

// unsigned long long 

template <> struct Number_type_traits<unsigned long long> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

inline
unsigned long long
div(unsigned long long i1, unsigned long long i2)
{ return i1 / i2; }

inline
double
to_double(unsigned long long i)
{ return static_cast<double>(i); }

inline
bool
is_finite(unsigned long long)
{ return true; }

inline
bool
is_valid(unsigned long long)
{ return true; }

namespace NTS {
  inline unsigned long long is_negative(unsigned long long i) 
    { return false; }
  
  inline unsigned long long sign(unsigned long long i) 
    { return is_positive(i) ? POSITIVE : ZERO; }
  
  inline unsigned long long is_positive(unsigned long long i)
    { return i != 0; }

  inline unsigned long long abs(unsigned long long i) { return i; }
} // namespace NTS

#endif // LONG_LONG

CGAL_END_NAMESPACE

#endif // CGAL_INT_H

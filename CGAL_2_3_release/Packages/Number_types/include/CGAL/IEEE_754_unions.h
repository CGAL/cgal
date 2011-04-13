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
// file          : IEEE_754_unions.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_IEEE_754_UNIONS_H
#define CGAL_IEEE_754_UNIONS_H

#include <iomanip>

union IEEE_754_double
{
  double   a;
#ifdef CGAL_BIG_ENDIAN
  struct { unsigned sign : 1;
           unsigned exp  :11;
           unsigned high :20;
           unsigned low  :32;
         } b;
  struct { unsigned H    :32;
           unsigned L    :32;
         } c;
#else
  struct { unsigned low  :32;
           unsigned sign : 1;
           unsigned exp  :11;
           unsigned high :20;
         } b;
  struct { unsigned L    :32;
           unsigned H    :32;
         } c;
#endif
};

union IEEE_754_float
{
  float    a;
  struct { unsigned sign : 1;
           unsigned exp  : 8;
           unsigned high :23;
         } b;
  unsigned c;
};

inline
void
show( IEEE_754_double* p)
{
  std::cout << std::endl;
  std::cout << std::hex << std::setw(8) << std::setfill('0') << p->c.H;
  std::cout << ' ';
  std::cout << std::hex << std::setw(8) << std::setfill('0') << p->c.L;
  std::cout << std::endl;
}

inline
void
show( IEEE_754_float* p)
{
  std::cout << std::endl;
  std::cout << std::hex << std::setw(8) << std::setfill('0') << p->c;
  std::cout << std::endl;
}

#endif // CGAL_IEEE_754_UNIONS_H

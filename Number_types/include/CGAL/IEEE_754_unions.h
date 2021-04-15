// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_IEEE_754_UNIONS_H
#define CGAL_IEEE_754_UNIONS_H

#include <iomanip>
#include <iostream>

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

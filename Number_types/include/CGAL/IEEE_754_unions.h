// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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

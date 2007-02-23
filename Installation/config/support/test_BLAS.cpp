// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
//
// Author(s)     : Laurent Saboret

// Test if BLAS is available

#include <iostream>
#include <cassert>

extern "C" {
  // Fortran interface...
  void dgemm (char*,char*,int*,
	      int*,int*,
	      double*,
	      double*,int*,double*,int*,
	      double*,
	      double*,int*);
  // ... or C interface
  void dgemm_(char*,char*,int*,
	      int*,int*,
	      double*,
	      double*,int*,double*,int*,
	      double*,
	      double*,int*);
}

namespace CGAL { namespace BLAS {

inline
void dgemm (char*a,char*b,int*c,
            int*d,int*e,
            double*f,
            double*g,int*h,double*i,int*j,
            double*k,
            double*l,int*m)
{
#ifdef CGAL_USE_F2C
  return ::dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m);
#else
  return ::dgemm (a,b,c,d,e,f,g,h,i,j,k,l,m);
#endif
}

} }

int main()
{
  char S[2] = "N";
  double A = 2.0;
  double B = 3.0;
  double C = 5.0;
  double alpha = 7.0;
  double beta  = 11.0;
  int n = 1;
  int ld = 1;

  CGAL::BLAS::dgemm(S,S, &n,&n,&n, &alpha, &A,&ld, &B,&ld, &beta, &C,&ld);
  assert(C == 97.0);

  std::cout << "ok for BLAS" << std::endl;

  return 0;
}

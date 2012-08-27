// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s)     : Marc Pouget

// Test if LAPACK is available

#include <iostream>
#include <cassert>

extern "C" {
  // taken from acml.h
void dgelss(int m, int n, int nrhs, 
            double *a, int lda, double *b, int ldb, double *sing, 
            double rcond, int *irank, int *info);

void dgelss_(int *m, int *n, int *nrhs,
                    double *a, int *lda, double *b, int *ldb, double *
                    s, double *rcond, int *rank, double *work, int *lwork,
                    int *info);
}

namespace CGAL { namespace LAPACK {

inline
void dgelss(int *m, int *n, int *nrhs,
       double *a, int *lda, double *b, int *ldb, double *
       s, double *rcond, int *rank, double *work, int *lwork,
       int *info)
{
#ifdef CGAL_USE_F2C
  ::dgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
#else
  ::dgelss(*m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank,  info);
#endif
}

} }

int main()
{
  double M=1, B=1;
  int m = 1,
    n = 1,
    nrhs = 1,
    lda = m,
    ldb = m,
    rank,
    lwork = 5*m,
    info;
  //c style
  double * sing_values = new double[n]; // * sing_values = (double*) malloc(n*sizeof(double));
  double* work = new double[lwork]; // (double*) malloc(lwork*sizeof(double));

  double rcond = -1;

  CGAL::LAPACK::dgelss(&m, &n, &nrhs, &M, &lda, &B, &ldb, sing_values, 
	  &rcond, &rank, work, &lwork, &info);
  assert(info==0);
  assert(B==1.);
  //clean up 
  delete sing_values;
  delete work;

  std::cout << "ok for lapack" << std::endl;
  
  return 0;
}

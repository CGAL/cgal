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
//
// Author(s)     : Laurent Saboret

// Test if BLAS is available

// blaswrap.h maps CBLAS function names to BLAS standard Fortran interface.
#ifdef CGAL_USE_CBLASWRAP
  #ifndef CGAL_USE_F2C
    #define CGAL_USE_F2C
  #endif
  #include <blaswrap.h>
#endif

#include <iostream>

extern "C" {
  // Fortran interface (taken from www.netlib.org/clapack)...
  void dgemm_(char* transa, char* transb, int* m,
              int* n, int* k,
              double* alpha,
              double* a, int* lda, double* b, int* ldb,
              double* beta,
              double* c, int* ldc);
  // ... or C interface (taken from AMD Core Math Library)
  void dgemm (char transa, char transb, int m,
              int n, int k,
              double alpha,
              double* a, int lda, double* b, int ldb,
              double beta,
              double* c, int ldc);
}

namespace CGAL { namespace BLAS {

inline
void dgemm (char* transa, char* transb, int* m,
            int* n, int* k,
            double* alpha,
            double* a, int* lda, double* b, int* ldb,
            double* beta,
            double* c, int* ldc)
{
#ifdef CGAL_USE_F2C
  ::dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  ::dgemm (*transa, *transb, *m, *n, *k, *alpha, a, *lda, b, *ldb, *beta, c, *ldc);
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
  if(C == 97.0) {
    std::cout << "ok for BLAS" << std::endl;
    return 0;
  }
  else {
    return 1;
  }
}

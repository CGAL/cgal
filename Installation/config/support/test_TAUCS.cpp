// Copyright (c) 2004  
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
// Author(s)     : Laurent Saboret

// Test if TAUCS is available


#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sys/types.h>
#include <sys/timeb.h>

#define TAUCS_CORE_DOUBLE

// In GCC 3.x, <complex.h> includes <complex> and
// complains if we include <complex.h> within "extern C {}"
#if defined(__GNUC__)
    #undef __DEPRECATED
    #include <complex.h>
#endif

extern "C"
{
    // Include TAUCS main header taucs.h
    #include <taucs.h>
}

// Avoid error with std::min() and std::max()
#undef min
#undef max


int main(int argc, char* argv[])
{
    // Create a TAUCS matrix to link with TAUCS main library
    int m = 4,n = 4,nnz = 4, i;
    taucs_ccs_matrix* pMatrix = taucs_ccs_create(m, n, nnz, TAUCS_DOUBLE);
    pMatrix->colptr[0] = 0;
    pMatrix->colptr[1] = 4;
    for ( i = 0; i < 4; i++ )
    {
        pMatrix->rowind[i] = i;
        pMatrix->taucs_values[i] = i;
    }

    // Call a method needing TAUCS external library Metis
    int*    perm;
    int*    invperm;
    taucs_ccs_order(pMatrix,
                    &perm,
                    &invperm,
                    (char*)"metis");

    // Call dpotrf() Lapack routine (Cholesky factorization)
    int sn_size = 0;
    int lda = 1;
    int info;
    taucs_potrf((char*)"LOWER",
                &sn_size,
                NULL,
                &lda,
                &info);

    // ftime() is in compat lib. on FreeBSD and crypto lib. on MacOSX
    struct timeb tb;
    ftime(&tb);

    // TAUCS provides no version number :-(
    // Version 1 is obsolete, thus we assume version 2 (latest is 2.2 on 03/2006)
    std::cout << "version=2.x" << std::endl;

    // clean up
    taucs_dccs_free(pMatrix);

    return 0;
}

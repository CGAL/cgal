// Copyright (c) 2004  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
//
// Author(s)     : Laurent Saboret

// ---------------------------------------------------------------------
// A short test program to evaluate a machine architecture.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

// Test if TAUCS is available
// First step "EXTERNTAUCS":
// - find TAUCS first include folder ${TAUCSROOT}/build/${TAUCS_OSTYPE}
// - find TAUCS external libs lapack f77blas cblas atlas metis g2c


#include <iostream>
#include <stdlib.h>
#include <stdio.h>

extern "C"
{
    // Include file in ${TAUCSROOT}/build/${TAUCS_OSTYPE}
    #include <taucs_config_build.h>

    /* from stuct.h in metis */
    typedef int idxtype;
    /* from metis.h */
    void METIS_NodeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *);
}


int main(int argc, char* argv[])
{
    // Call function in metis lib
    METIS_NodeND(NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    // TAUCS provides no version number :-(
    // Version 1 is obsolete, thus we assume version 2 (latest is 2.2 on 03/2006)
    std::cout << "version=2" << std::endl;

    return 0;
}



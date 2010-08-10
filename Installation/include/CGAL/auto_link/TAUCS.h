// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
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

#ifndef CGAL_AUTO_LINK_TAUCS_H
#define CGAL_AUTO_LINK_TAUCS_H

// Skip the whole file if auto-link is off
#if !defined(CGAL_NO_AUTOLINK_TAUCS) && !defined(CGAL_NO_AUTOLINK)


#if defined(_WIN32) && !defined(_WIN64) // if Windows 32 bits

// CGAL ships with TAUCS for Windows 32 bits, i.e this set of libraries (e.g. for VC++ 8 /MD):
// libtaucs-vc80-mt.lib libmetis-vc80-mt.lib liblapack.lib libf77blas.lib libcblas.lib libatlas.lib vcf2c-vc80-mt.lib
//
// Notes: - Order matters.
//        - Libraries with no "vc" toolset are compiled by gcc/g77. They are
//          compatible with VC++ 7.1, 8.0 and 9.0, and with all VC++ runtimes.
//        - Tested with 7.1, 8.0 and 9.0.

#define CGAL_LIB_NAME libtaucs
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME libmetis
#include <CGAL/auto_link/auto_link.h>

// Link with LAPACK, BLAS and F2C
#include <CGAL/auto_link/LAPACK.h>

#endif // Win32


#ifdef _WIN64 // if Windows 64 bits

// ATLAS is not compatible with Win64, thus CGAL ships with CLAPACK and CBLAS.
// VC++ >= 8.0 is compatible with Windows 64 bits.
// The set set of libraries is (e.g. for VC++ 8 /MD):
// libtaucs-vc80-mt.lib libmetis-vc80-mt.lib clapack-vc80-mt.lib blas-vc80-mt.lib vcf2c-vc80-mt.lib
//
// Notes: - Order matters.
//        - Tested with VC++ 8.0 and 9.0.

#define CGAL_LIB_NAME libtaucs
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME libmetis
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

// Link with LAPACK, BLAS and F2C
#include <CGAL/auto_link/LAPACK.h>

#endif // _WIN64


#endif // CGAL_NO_AUTOLINK_TAUCS && CGAL_NO_AUTOLINK

#endif // CGAL_AUTO_LINK_TAUCS_H

// Copyright (c) 1997-2004  
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Laurent Saboret

#ifndef CGAL_AUTO_LINK_LAPACK_H
#define CGAL_AUTO_LINK_LAPACK_H

#include <CGAL/config.h>

// Skip the whole file if auto-link is off
#if !defined(CGAL_NO_AUTOLINK_LAPACK) && !defined(CGAL_NO_AUTOLINK)


#if defined(_WIN32) && !defined(_WIN64) // if Windows 32 bits

// CGAL ships with ATLAS for Windows 32 bits, i.e this set of libraries (e.g. for VC++ 8 /MD):
// liblapack.lib libf77blas.lib libcblas.lib libatlas.lib vcf2c-vc80-mt.lib
//
// Notes: - Order matters.
//        - Libraries with no "vc" toolset are compiled by gcc/g77. They are
//          compatible with VC++ 7.1, 8.0 and 9.0, and with all VC++ runtimes.
//        - Tested with 7.1, 8.0 and 9.0.

#define CGAL_LIB_NAME liblapack
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME libf77blas
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME libcblas
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME libatlas
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME vcf2c
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

// ATLAS provides BLAS and LAPACK standard Fortran interface
#ifndef CGAL_USE_F2C
  #define CGAL_USE_F2C
#endif

#endif // Win32


#ifdef _WIN64 // if Windows 64 bits

// ATLAS is not compatible with Win64, thus CGAL ships with CLAPACK and CBLAS.
// VC++ >= 8.0 is compatible with Windows 64 bits.
// The set set of libraries is (e.g. for VC++ 8 /MD):
// clapack-vc80-mt.lib blas-vc80-mt.lib vcf2c-vc80-mt.lib
//
// Notes: - Order matters.
//        - Tested with VC++ 8.0 and 9.0.

#define CGAL_LIB_NAME clapack
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME blas
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

#define CGAL_LIB_NAME vcf2c
#define CGAL_AUTO_LINK_NOMANGLE
#include <CGAL/auto_link/auto_link.h>

// CLAPACK provides LAPACK standard Fortran interface.
// blaswrap.h maps CBLAS function names to BLAS standard Fortran interface.
#ifndef CGAL_USE_F2C
  #define CGAL_USE_F2C
#endif
#include <blaswrap.h>

#endif // _WIN64


#endif // CGAL_NO_AUTOLINK_LAPACK && CGAL_NO_AUTOLINK

#endif // CGAL_AUTO_LINK_LAPACK_H

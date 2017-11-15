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

#ifndef CGAL_AUTO_LINK_TAUCS_H
#define CGAL_AUTO_LINK_TAUCS_H

#include <CGAL/config.h>

// Skip the whole file if auto-link is off
#if !defined(CGAL_NO_AUTOLINK_TAUCS) && !defined(CGAL_NO_AUTOLINK)

#  if defined(_WIN32) || defined(_WIN64) 

#    define CGAL_LIB_NAME libtaucs
#    define CGAL_AUTO_LINK_NOMANGLE
#    include <CGAL/auto_link/auto_link.h>

#    define CGAL_LIB_NAME libmetis
#    define CGAL_AUTO_LINK_NOMANGLE
#    include <CGAL/auto_link/auto_link.h>

// Link with LAPACK, BLAS and F2C
#    include <CGAL/auto_link/LAPACK.h>

#  endif // Win32|Win64

#endif // CGAL_NO_AUTOLINK_TAUCS && CGAL_NO_AUTOLINK

#endif // CGAL_AUTO_LINK_TAUCS_H

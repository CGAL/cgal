// Copyright (c) 2016  GeometryFactory (France). All rights reserved.
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
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H
#define CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H

// GMP and MPFR are highly recommended in CGAL.
#define CGAL_USE_GMP 1
#define CGAL_USE_MPFR 1

#if CGAL_DISABLE_GMP && ! defined(CGAL_NO_GMP)
#  define CGAL_NO_GMP 1
#endif

#if CGAL_NO_GMP || CGAL_NO_MPFR
#  undef CGAL_USE_MPFR
#  undef CGAL_USE_GMP
#endif

#if defined(__has_include)
#  if CGAL_USE_GMP && ! __has_include(<gmp.h>)
#    warning "<gmp.h> cannot be found. Less efficient number types will be used instead. Define CGAL_NO_GMP=1 if that is on purpose."
#    undef CGAL_USE_GMP
#    undef CGAL_USE_MPFR
#  elif CGAL_USE_MPFR && ! __has_include(<mpfr.h>)
#    warning "<mpfr.h> cannot be found and the GMP support in CGAL requires it. Less efficient number types will be used instead. Define CGAL_NO_GMP=1 if that is on purpose."
#    undef CGAL_USE_GMP
#    undef CGAL_USE_MPFR
#  endif // CGAL_USE_MPFR and no <mpfr.h>
#endif // __has_include

#if CGAL_USE_GMP && CGAL_USE_MPFR && ! CGAL_NO_CORE
#  define CGAL_USE_CORE 1
#endif

#endif // CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H

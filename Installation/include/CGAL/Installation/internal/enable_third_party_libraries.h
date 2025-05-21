// Copyright (c) 2016  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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

#if defined(__has_include) && ( ! defined _MSC_VER || _MSC_VER > 1900)
#  if CGAL_USE_GMP && ! __has_include(<gmp.h>)
#    undef CGAL_USE_GMP
#    undef CGAL_USE_MPFR
#  elif CGAL_USE_MPFR && ! __has_include(<mpfr.h>)
#    undef CGAL_USE_GMP
#    undef CGAL_USE_MPFR
#  endif // CGAL_USE_MPFR and no <mpfr.h>
#endif // __has_include


// It is easier to disable this number type completely for old versions.
// Before 1.63, I/O is broken.  Again, disabling the whole file is just the
// easy solution.
// MSVC had trouble with versions <= 1.69:
// https://github.com/boostorg/multiprecision/issues/98
//
// Disable also on Windows 32 bits
// because CGAL/cpp_float.h assumes _BitScanForward64 is available
// See https://learn.microsoft.com/en-us/cpp/intrinsics/bitscanforward-bitscanforward64
//
// Disable also with PowerPC processors, with Boost<1.80 because of that bug:
// https://github.com/boostorg/multiprecision/pull/421
//
#if !defined CGAL_DO_NOT_USE_BOOST_MP && \
    (!defined _MSC_VER || BOOST_VERSION >= 107000) && \
    (!defined _WIN32 || defined _WIN64) && \
    (BOOST_VERSION >= 108000 || (!defined _ARCH_PPC && !defined _ARCH_PPC64))
#define CGAL_USE_BOOST_MP 1
#endif


#if CGAL_USE_BOOST_MP
#if ! CGAL_NO_CORE
#  define CGAL_USE_CORE 1
#endif
#endif

#endif // CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H

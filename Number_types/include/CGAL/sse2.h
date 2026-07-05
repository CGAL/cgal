// Copyright (c) 1999,2007
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Andreas Fabri

#ifndef CGAL_SSE2_H
#define CGAL_SSE2_H


#if defined ( _MSC_VER )
#define CGAL_ALIGN_16  __declspec(align(16))
#elif defined( __GNUC__ ) || defined(__clang__)
#define  CGAL_ALIGN_16 __attribute__((aligned(16)))
#endif

// Architecture-specific SIMD headers
#if defined(__aarch64__) || defined(_M_ARM64)
#  include <arm_neon.h>
#else
#  include <emmintrin.h>
#endif

#endif // CGAL_SSE2_H

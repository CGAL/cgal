// Copyright (c) 2010-2011  GeometryFactory Sarl (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau


// The main goal of FPU.h is to define functions or macros to modify the
// control word of the FPU, to:
//   - set the precision to 53 bits of mantissa,
//   - get/set the rounding mode.
//
// The goal of FPU_extension.h is to define inline functions similar to
// feclearexcept and fetestexcept of C99.
//
// For the moment, only i386 and x64 processors are supported, with MSVC,
// gcc, or the Intel compiler suite. Otherwise, the non-inline functions of
// C99 are used.

#ifndef CGAL_FPU_EXTENSION_H
#define CGAL_FPU_EXTENSION_H

#if __i386__ && !defined __PGI && !defined __SUNPRO_CC
#  ifdef CGAL_SAFE_SSE2
#    include <CGAL/FPU_gcc_i386_sse2.h>
#  else
#    include <CGAL/FPU_gcc_i386.h>
#  endif
#elif defined _MSC_VER
#  include <CGAL/FPU_msvc.h>
#else

// generic functions, using C99

extern "C" {
#  include <fenv.h>
}

namespace CGAL {

inline int
feclearexcept(int exceptions) {
  return ::feclearexcept(exceptions);
}

inline int
fetestexcept(int exceptions) {
  return ::fetestexcept(exceptions);
}

} // end namespace CGAL

#endif // use fenv

#endif // CGAL_FPU_EXTENSION_H

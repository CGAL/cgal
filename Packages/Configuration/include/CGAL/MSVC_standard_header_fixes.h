// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/MSVC_standard_header_fixes.h
// chapter       : $CGAL_Chapter: Configuration $
//
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Utrecht University
// ============================================================================

#ifndef CGAL_MSVC_STANDARD_HEADER_FIXES_H
#define CGAL_MSVC_STANDARD_HEADER_FIXES_H

#pragma warning(once: 4291)
#pragma warning(once:4503)


// the iterator specializations should be used for:
// cl 1300 and Intel Compiler
#if defined (_MSC_VER) && ( (_MSC_VER == 1300) || defined (__INTEL_COMPILER) )
#include <CGAL/config/msvc7/iterator_specializations.h>
#endif

#include <cmath>
namespace std {
	using ::sqrt;
}

#include <cstddef>
namespace std{
  using ::size_t;
  using ::ptrdiff_t;
}

#include <ctime>
namespace std{
  using ::time_t;
}

#endif // CGAL_MSVC_STANDARD_HEADER_FIXES_H

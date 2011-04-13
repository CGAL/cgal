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

#include <cmath>
namespace std {
	using ::sqrt;
}
#define M_PI 3.14159265358979323846

#include <cstddef>
namespace std{
  using ::size_t;
  using ::ptrdiff_t;
}

#endif // CGAL_MSVC_STANDARD_HEADER_FIXES_H

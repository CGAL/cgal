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

#include <cmath>
namespace std {
	using ::fabs;
	using ::sqrt;
	using ::log;
	using ::cos;
	using ::sin;
	using ::ceil;
	using ::floor;
}
#define M_PI 3.14159
#include <cstddef>
namespace std{
using ::size_t;
using ::ptrdiff_t;
}

#include <cstdlib>
namespace std{
using ::abort;
using ::atoi;
}
#include <cstring>
namespace std{
using ::strcat;
using ::strcpy;
}

#include <ctime>

namespace std{
using ::clock;
using ::clock_t;
}

#endif // CGAL_MSVC_STANDARD_HEADER_FIXES_H

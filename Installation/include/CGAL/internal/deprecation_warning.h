// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Philipp Möller, Mael Rouxel-Labbé

// Including this header will cause compilation to fail
// if CGAL_NO_DEPRECATED_CODE is defined. If this is not the case, it will issue
// a warning during compilation, unless CGAL_NO_DEPRECATION_WARNINGS is defined.

// CGAL_DEPRECATED_HEADER, CGAL_REPLACEMENT_HEADER, and
// CGAL_DEPRECATED_MESSAGE_DETAILS can be defined
// to a string literal to customize the warning.

// CGAL_DEPRECATED_HEADER, CGAL_REPLACEMENT_HEADER, and
// CGAL_DEPRECATED_MESSAGE_DETAILS are undefined after the file is included.

// The lack of an include guard is intentional and necessary.

#include <CGAL/assertions.h>

// whether to print Warning or Error
#if defined(CGAL_NO_DEPRECATED_CODE)
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE_STATUS "Error: "
#  define CGAL_INTERNAL_NO_DEPRECATED_CODE_MESSAGE " and CGAL_NO_DEPRECATED_CODE is defined."
#else
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE_STATUS "Warning: "
#  define CGAL_INTERNAL_NO_DEPRECATED_CODE_MESSAGE "."
#endif

// if the name of the deprecated header is given, print it
#if defined(CGAL_DEPRECATED_HEADER)
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE_DEPRECATED_HEADER \
     "The header `" CGAL_DEPRECATED_HEADER "` is deprecated"
#else
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE_DEPRECATED_HEADER \
     "A deprecated header has been included"
#endif

// if a replacement header is given, print it
#if defined(CGAL_REPLACEMENT_HEADER)
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE_HEADERS \
            CGAL_INTERNAL_DEPRECATED_MESSAGE_DEPRECATED_HEADER \
            CGAL_INTERNAL_NO_DEPRECATED_CODE_MESSAGE \
            " Please use `" CGAL_REPLACEMENT_HEADER "` instead. "
#else
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE_HEADERS \
            CGAL_INTERNAL_DEPRECATED_MESSAGE_DEPRECATED_HEADER \
            CGAL_INTERNAL_NO_DEPRECATED_CODE_MESSAGE " "
#endif

// if more details are given, print them
#if defined(CGAL_DEPRECATED_MESSAGE_DETAILS)
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE \
     CGAL_INTERNAL_DEPRECATED_MESSAGE_STATUS \
     CGAL_INTERNAL_DEPRECATED_MESSAGE_HEADERS \
     "Additional information: "\
     CGAL_DEPRECATED_MESSAGE_DETAILS
#else
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE \
     CGAL_INTERNAL_DEPRECATED_MESSAGE_STATUS \
     CGAL_INTERNAL_DEPRECATED_MESSAGE_HEADERS
#endif

#if defined(CGAL_NO_DEPRECATED_CODE) // No deprecated code.
CGAL_static_assertion_msg(false, CGAL_INTERNAL_DEPRECATED_MESSAGE);
#elif !defined(CGAL_NO_DEPRECATION_WARNINGS) // don't trigger on NO_DEPRECATION_WARNINGS
#  if defined(_MSC_VER) || defined(__BORLANDC__) || defined(__DMC__)
#    pragma message (CGAL_INTERNAL_DEPRECATED_MESSAGE)
#  elif defined(__GNUC__) || defined(__HP_aCC) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
     // warning does not expand its arguments, issue a warning and add the message.
#    warning "A deprecated header has been included."
#    pragma message (CGAL_INTERNAL_DEPRECATED_MESSAGE)
#  endif //defined
#endif

// those macros have been defined in all cases
#undef CGAL_INTERNAL_DEPRECATED_MESSAGE_STATUS
#undef CGAL_INTERNAL_DEPRECATED_MESSAGE_DEPRECATED_HEADER
#undef CGAL_INTERNAL_DEPRECATED_MESSAGE_HEADERS
#undef CGAL_INTERNAL_DEPRECATED_MESSAGE

#if defined(CGAL_DEPRECATED_MESSAGE_DETAILS)
#  undef CGAL_DEPRECATED_MESSAGE_DETAILS
#endif

#if defined(CGAL_DEPRECATED_HEADER)
#  undef CGAL_DEPRECATED_HEADER
#endif

#if defined(CGAL_REPLACEMENT_HEADER)
#  undef CGAL_REPLACEMENT_HEADER
#endif

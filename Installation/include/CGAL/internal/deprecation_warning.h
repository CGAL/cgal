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
// 
//
// Author(s)     : Philipp Möller

// Including this header will issue a warning during compilation or
// cause compilation to fail if CGAL_NO_DEPRECATED_CODE is defined.
// CGAL_DEPRECATED_HEADER and CGAL_REPLACEMENT_HEADER can be defined
// to a string literal to customize the warning.
// CGAL_DEPRECATED_HEADER and CGAL_REPLACEMENT_HEADER are undefined,
// after the file is included.  

// The lack of an include guard is intentional and necessary.

#include <CGAL/assertions.h>

#ifndef CGAL_NO_DEPRECATION_WARNINGS

#if defined(CGAL_NO_DEPRECATED_CODE)
// No deprecated code.
CGAL_static_assertion_msg(false, "A deprecated header has been included and CGAL_NO_DEPRECATED_CODE is defined.");
#endif // CGAL_NO_DEPRECATED_CODE

// Build the message
#define CGAL_INTERNAL_DEPRECATED_MESSAGE "Warning: A deprecated header has been included."

#if defined(CGAL_DEPRECATED_HEADER) && defined(CGAL_REPLACEMENT_HEADER)
#  undef CGAL_INTERNAL_DEPRECATED_MESSAGE
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE "Warning: The header " CGAL_DEPRECATED_HEADER " is deprecated. " \
                                           "Please use " CGAL_REPLACEMENT_HEADER " instead."
#elif defined(CGAL_DEPRECATED_HEADER)
#  undef CGAL_INTERNAL_DEPRECATED_MESSAGE
#  define CGAL_INTERNAL_DEPRECATED_MESSAGE "Warning: The header " CGAL_DEPRECATED_HEADER " is deprecated."
#endif

// don't trigger on NO_DEPRECATIOON_WARNINGS and don't trigger twice on NO_DEPRECATED_CODE
#if !defined(CGAL_NO_DEPRECATION_WARNINGS) && !defined(CGAL_NO_DEPRECATED_CODE)
#  if defined(_MSC_VER) || defined(__BORLANDC__) || defined(__DMC__)
#    pragma message (CGAL_INTERNAL_DEPRECATED_MESSAGE)
#  elif defined(__GNUC__) || defined(__HP_aCC) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
     // warning does not expand its arguments, issue a warning and add the message.
#    warning "A deprecated header has been included."
#    pragma message (CGAL_INTERNAL_DEPRECATED_MESSAGE)
#  endif //defined
#endif

#endif // CGAL_NO_DEPRECATION_WARNINGS

#if defined(CGAL_DEPRECATED_HEADER)
#  undef CGAL_DEPRECATED_HEADER
#endif

#if defined(CGAL_REPLACEMENT_HEADER)
#  undef CGAL_REPLACEMENT_HEADER
#endif

#if defined(CGAL_INTERNAL_DEPRECATED_MESSAGE)
#  undef CGAL_INTERNAL_DEPRECATED_MESSAGE
#endif

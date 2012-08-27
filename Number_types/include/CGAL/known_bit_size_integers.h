// Copyright (c) 1999-2005  
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
//
//
// Author(s)     : Lutz Kettner, Sylvain Pion

// This header has been deprecated and all types are aliases for the
// respective type in boost/cstdint.hpp.
// boost/cstdint.hpp is preferred over cstdint as there are numerous
// compilers and platforms with slight variations in availability of
// the C99 headers and their C++0x equivalents

#ifndef CGAL_KNOWN_BIT_SIZE_INTEGERS_H
#define CGAL_KNOWN_BIT_SIZE_INTEGERS_H

#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(__DMC__)
#  pragma message ("Warning: This header is deprecated. Please use boost/cstdint.hpp instead")
#elif defined(__GNUC__) || defined(__HP_aCC) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
#  warning "This header is deprecated. Please use boost/cstdint.hpp instead"
#endif

#include <CGAL/number_type_basic.h>
#include <boost/cstdint.hpp>

#ifndef CGAL_NO_DEPRECATED_CODE

namespace CGAL {

  typedef boost::int8_t    Integer8;
  typedef boost::int16_t   Integer16;
  typedef boost::int32_t   Integer32;

  typedef boost::uint8_t   UInteger8;
  typedef boost::uint16_t  UInteger16;
  typedef boost::uint32_t  UInteger32;

#ifndef BOOST_NO_INT64_T
// this macro is still provided but its use is discouraged
#   define CGAL_HAS_INTEGER64
  typedef boost::int64_t   Integer64;
  typedef boost::uint64_t  UInteger64;
#endif

} //namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_KNOWN_BIT_SIZE_INTEGERS_H

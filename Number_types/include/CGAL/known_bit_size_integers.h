// Copyright (c) 1999-2005
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
// Author(s)     : Lutz Kettner, Sylvain Pion

// This header has been deprecated and all types are aliases for the
// respective type in boost/cstdint.hpp.
// boost/cstdint.hpp is preferred over cstdint as there are numerous
// compilers and platforms with slight variations in availability of
// the C99 headers and their C++0x equivalents

#ifndef CGAL_KNOWN_BIT_SIZE_INTEGERS_H
#define CGAL_KNOWN_BIT_SIZE_INTEGERS_H

#define CGAL_DEPRECATED_HEADER "<CGAL/known_bit_size_integers.h>"
#define CGAL_REPLACEMENT_HEADER "<boost/cstdint.hpp>"
#include <CGAL/internal/deprecation_warning.h>

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

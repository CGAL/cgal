// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_DEMANGLE_H
#define CGAL_DEMANGLE_H

#if BOOST_VERSION >= 105600
#include <boost/core/demangle.hpp>
#else
#include <boost/units/detail/utility.hpp>
#endif

namespace CGAL {


inline std::string demangle(const char* name)
{
#if BOOST_VERSION >= 105600
  return boost::core::demangle(name);
#else
  return boost::units::detail::demangle(name);
#endif
}


} //namespace CGAL

#endif // CGAL_DEMANGLE_H

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

#include <boost/core/demangle.hpp>

namespace CGAL {


inline std::string demangle(const char* name)
{
  return boost::core::demangle(name);
}


} //namespace CGAL

#endif // CGAL_DEMANGLE_H

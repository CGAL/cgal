// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_FOREACH_H
#define CGAL_FOREACH_H

#include <CGAL/config.h>
#include <CGAL/Iterator_range.h>

#ifndef DOXYGEN_RUNNING
#if BOOST_VERSION >= 105100 && !defined(BOOST_NO_CXX11_RANGE_BASED_FOR)

namespace CGAL{
namespace for_range_loop{

template <class T>
Iterator_range<T>
constexpr make_range(const std::pair<T,T>& p)
{
  return CGAL::make_range(p);
}

template <class T>
Iterator_range<T>
constexpr make_range(std::pair<T,T>& p)
{
  return CGAL::make_range(p);
}

template <class T>
Iterator_range<T>
constexpr make_range(std::pair<T,T>&& p)
{
  return CGAL::make_range(p);
}

template <class T>
constexpr T& make_range(T& t)
{
  return t;
}

template <class T>
constexpr T make_range(T&& t)
{
  return std::move(t);
}

} }

#define CGAL_FOREACH(A,B) for(A : CGAL::for_range_loop::make_range(B))
#else
#include <boost/foreach.hpp>
#define CGAL_FOREACH(A,B) BOOST_FOREACH(A, B)
#endif
#else
/// \ingroup  PkgStlExtension
/// If the version of `boost` is at least 1.51 and the compiler used support range-based for loops
/// introduced with \cpp11, use it; otherwise fallback onto `BOOST_FOREACH`.
#define CGAL_FOREACH(A,B)
#endif
#endif

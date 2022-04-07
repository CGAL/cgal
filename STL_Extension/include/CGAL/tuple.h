// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot, Sylvain Pion

#ifndef CGAL_TUPLE_H
#define CGAL_TUPLE_H

#  include <tuple>

namespace CGAL {

// Tool to test whether a type V is among the types of a tuple<...> = T.
template <typename V, typename T>
struct Is_in_tuple;

template <typename V, typename T0, typename... T>
struct Is_in_tuple <V, std::tuple<T0, T...> >
{
  static const bool value = Is_in_tuple<V, std::tuple<T...> >::value;
};

template <typename V, typename... T>
struct Is_in_tuple <V, std::tuple<V, T...> >
{
  static const bool value = true;
};

template <typename V>
struct Is_in_tuple <V, std::tuple<> >
{
  static const bool value = false;
};

} //namespace CGAL

#endif // CGAL_TUPLE_H

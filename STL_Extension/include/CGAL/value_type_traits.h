// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Alberto Ganesh Barbati and Laurent Saboret

#ifndef CGAL_VALUE_TYPE_TRAITS_H
#define CGAL_VALUE_TYPE_TRAITS_H

#include <iterator>

namespace CGAL {

/// \ingroup  PkgSTLExtensionRef
/// Class providing the value type of an iterator, and
/// in the case of an output iterator, a type of objects that can be put in it.
///
template <class T>
struct value_type_traits
{
  #ifndef DOXYGEN_RUNNING
  typedef typename std::iterator_traits<T>::value_type type;
  #else
  /// If `T` is `std::insert_iterator<Container>`, `std::back_insert_iterator<Container>` or
  /// `std::front_insert_iterator<Container>`, then `type` is `Container::value_type`.
  /// Otherwise, `type` is `std::iterator_traits<T>::%value_type`.

  typedef unspecified_type type;
  #endif
};

template <class Container>
struct value_type_traits<std::back_insert_iterator<Container> >
{
  typedef typename Container::value_type type;
};

template <class Container>
struct value_type_traits<std::insert_iterator<Container> >
{
  typedef typename Container::value_type type;
};

template <class Container>
struct value_type_traits<std::front_insert_iterator<Container> >
{
  typedef typename Container::value_type type;
};

} //namespace CGAL

#endif // CGAL_VALUE_TYPE_TRAITS_H

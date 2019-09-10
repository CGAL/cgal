// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
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
// Author(s) : Alberto Ganesh Barbati and Laurent Saboret

#ifndef CGAL_VALUE_TYPE_TRAITS_H
#define CGAL_VALUE_TYPE_TRAITS_H

#include <iterator>

namespace CGAL {

/// \ingroup  PkgStlExtension
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

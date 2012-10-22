// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s) : Alberto Ganesh Barbati and Laurent Saboret

#ifndef CGAL_VALUE_TYPE_TRAITS_H
#define CGAL_VALUE_TYPE_TRAITS_H

#include <iterator>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

/// Traits class to get the value type of any iterator,
/// including an output iterator.
/// Based on code posted by Alberto Ganesh Barbati at
/// http://www.adras.com/Why-no-std-back-insert-iterator-value-type.t2639-153-3.html
///
/// Usage is:
/// typedef typename value_type_traits<Iter>::type value_type;

template <class T>
struct value_type_traits
{
  typedef typename std::iterator_traits<T>::value_type type;
};

template <class T>
struct value_type_traits<std::back_insert_iterator<T> >
{
  typedef typename T::value_type type;
};

/// \endcond

} //namespace CGAL

#endif // CGAL_VALUE_TYPE_TRAITS_H

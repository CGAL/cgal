// Copyright (c) 2008-2009 GeometryFactory and INRIA
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Andreas Fabri and laurent Saboret

#ifndef CGAL_POINT_SET_PROPERTY_MAP_H
#define CGAL_POINT_SET_PROPERTY_MAP_H

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/value_type_traits.h>

#include <boost/property_map.hpp>
#include <boost/tuple/tuple.hpp>

#include <utility> // defines std::pair

namespace CGAL {


  //=========================================================================
  /// Property map T* -> T.
  /// A common usage is a property map Point_with_normal_3* -> position (Point_3).

  template <typename T>
  struct Dereference_property_map
    : public boost::put_get_helper<T&, Dereference_property_map<T> >
  {
    typedef T* key_type;
    typedef T value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    /// Access the map elements.
    template <class Iter> // Type convertible to key_type
    reference operator[](Iter it) const { return reference(*it); }
  };

  /// Free function to create a Dereference_property_map property map
  template <class Iter> // Type convertible to key_type
  Dereference_property_map<typename CGAL::value_type_traits<Iter>::type>
  make_dereference_property_map(Iter)
  {
    // value_type_traits is a workaround as back_insert_iterator's value_type is void
    return Dereference_property_map<typename CGAL::value_type_traits<Iter>::type>();
  }


  //=========================================================================
  /// Property map Point_with_normal_3* -> normal vector (Vector_3)

  template <class Gt>
  struct Normal_vector_property_map
    : public boost::put_get_helper<typename Gt::Vector_3&,
                                   Normal_vector_property_map<Gt> >
  {
    typedef Point_with_normal_3<Gt> Point_with_normal; ///< Position + normal
    typedef typename Gt::Vector_3 Vector; /// normal

    typedef Point_with_normal* key_type;
    typedef Vector value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    /// Access the map elements.
    template <class Iter> // Type convertible to key_type
    reference operator[](Iter it) const { return (reference) it->normal(); }
  };

  /// Free function to create a Normal_vector_property_map property map
  template <class Iter> // Type convertible to key_type
  Normal_vector_property_map<typename CGAL::Kernel_traits<typename CGAL::value_type_traits<Iter>::type>::Kernel>
  make_normal_vector_property_map(Iter)
  {
    // value_type_traits is a workaround as back_insert_iterator's value_type is void
    typedef typename CGAL::value_type_traits<Iter>::type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    return Normal_vector_property_map<Kernel>();
  }


  //=========================================================================
  // Property maps Pair* -> Pair::first_type
  // and           Pair* -> Pair::second_type.

  /// Property map Pair* -> Pair::first_type
  template <typename Pair>
  struct First_of_pair_property_map
    : public boost::put_get_helper<typename Pair::first_type&,
                                   First_of_pair_property_map<Pair> >
  {
    typedef Pair* key_type;
    typedef typename Pair::first_type value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    /// Access the map elements.
    template <class Iter> // Type convertible to key_type
    reference operator[](Iter pair) const { return reference(pair->first); }
  };

  /// Free function to create a First_of_pair_property_map property map
  template <class Iter> // Type convertible to key_type
  First_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>
  make_first_of_pair_property_map(Iter)
  {
    // value_type_traits is a workaround as back_insert_iterator's value_type is void
    return First_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>();
  }

  /// Property map Pair* -> Pair::second_type
  template <typename Pair>
  struct Second_of_pair_property_map
    : public boost::put_get_helper<typename Pair::second_type&,
                                   Second_of_pair_property_map<Pair> >
  {
    typedef Pair* key_type;
    typedef typename Pair::second_type value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    /// Access the map elements.
    template <class Iter> // Type convertible to key_type
    reference operator[](Iter pair) const { return reference(pair->second); }
  };

  /// Free function to create a Second_of_pair_property_map property map
  template <class Iter> // Type convertible to key_type
  Second_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>
  make_second_of_pair_property_map(Iter)
  {
    // value_type_traits is a workaround as back_insert_iterator's value_type is void
    return Second_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>();
  }


  //=========================================================================
  /// Property map Tuple* -> Nth element of tuple

  template <int N, typename Tuple>
  struct Nth_of_tuple_property_map
    : public boost::put_get_helper<typename boost::tuples::element<N,Tuple>::type&,
                                   Nth_of_tuple_property_map<N,Tuple> >
  {
    typedef Tuple* key_type;
    typedef typename boost::tuples::element<N,Tuple>::type value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    /// Access the map elements.
    template <class Iter> // Type convertible to key_type
    reference operator[](Iter tuple) const { return (reference) tuple->get<N>(); }
  };

  /// Free function to create a Nth_of_tuple_property_map property map
  template <int N, class Iter> // Type convertible to key_type
  Nth_of_tuple_property_map<N, typename CGAL::value_type_traits<Iter>::type>
  make_nth_of_tuple_property_map(Iter)
  {
    // value_type_traits is a workaround as back_insert_iterator's value_type is void
    return Nth_of_tuple_property_map<N, typename CGAL::value_type_traits<Iter>::type>();
  }


} // namespace CGAL

#endif // CGAL_POINT_SET_PROPERTY_MAP_H

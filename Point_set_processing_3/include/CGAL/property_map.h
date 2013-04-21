// Copyright (c) 2008-2009 GeometryFactory and INRIA
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
// Author(s)     : Andreas Fabri and Laurent Saboret

#ifndef CGAL_POINT_SET_PROPERTY_MAP_H
#define CGAL_POINT_SET_PROPERTY_MAP_H

#include <CGAL/value_type_traits.h>

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
  #include <boost/property_map/property_map.hpp>
#else
  #include <boost/property_map.hpp>
#endif
#include <boost/tuple/tuple.hpp>

#include <utility> // defines std::pair

namespace CGAL {

/// An alternative to boost::put_get_helper for:
/// - return type of `get` function is `const Reference`
/// - `key` is passed by reference in `put` function
template <class Reference, class LvaluePropertyMap>
struct put_get_helper_pass_key_by_reference { };

template <class PropertyMap, class Reference, class K>
inline const Reference
get(const put_get_helper_pass_key_by_reference<Reference, PropertyMap>& pa, const K& k)
{
  Reference v = static_cast<const PropertyMap&>(pa)[k];
  return v;
}

template <class PropertyMap, class Reference, class K, class V>
inline void
put(const put_get_helper_pass_key_by_reference<Reference, PropertyMap>& pa, K& k, const V& v)
{
  static_cast<const PropertyMap&>(pa)[k] = v;
}


/// \ingroup PkgProperty_map
/// Property map that converts a `T*` pointer (or in general an iterator
/// over `T` elements) to the `T` object.
///
/// \cgalModels `LvaluePropertyMap`
template <typename T>
struct Dereference_property_map
  : public boost::put_get_helper<T&, Dereference_property_map<T> >
{
  typedef T* key_type; ///< typedef to 'T*'
  typedef T value_type; ///< typedef to 'T'
  typedef value_type& reference; ///< typedef to 'T&'
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`

  /// Access a property map element.
  ///
  /// @tparam Iter Type convertible to `key_type`.
  template <class Iter>
  reference operator[](Iter it) const { return reference(*it); }
};

/// Free function to create a `Dereference_property_map` property map.
///
/// \relates Dereference_property_map 
template <class Iter> // Type convertible to `key_type`
Dereference_property_map<typename CGAL::value_type_traits<Iter>::type>
make_dereference_property_map(Iter)
{
  // value_type_traits is a workaround as back_insert_iterator's `value_type` is void
  return Dereference_property_map<typename CGAL::value_type_traits<Iter>::type>();
}

template <typename T>
struct Typed_identity_property_map_by_reference
  : put_get_helper_pass_key_by_reference<T, Typed_identity_property_map_by_reference<T> >
{
  typedef T key_type;
  typedef T value_type;
  typedef T& reference;
  typedef boost::lvalue_property_map_tag category;

  reference operator[](key_type& v) const { return reference(v); }
  const reference operator[](const key_type& v) const { return reference(v); }
};



//=========================================================================
#ifdef CGAL_USE_OLD_PAIR_PROPERTY_MAPS
/// \ingroup PkgProperty_map
/// Property map that accesses the first item of a `std::pair`. 
/// \tparam Pair Instance of `std::pair`. 
/// \cgalModels `LvaluePropertyMap`
///
/// \sa `CGAL::Second_of_pair_property_map<Pair>`
template <typename Pair>
struct First_of_pair_property_map
  : public boost::put_get_helper<typename Pair::first_type&,
                                 First_of_pair_property_map<Pair> >
{
  typedef Pair* key_type; ///< typedef to 'Pair*'
  typedef typename Pair::first_type value_type; ///< typedef to `Pair::first_type`
  typedef value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< boost::lvalue_property_map_tag

  /// Access a property map element.
  ///
  /// @tparam Iter Type convertible to `key_type`.
  template <class Iter>
  reference operator[](Iter pair) const { return reference(pair->first); }
};
#else
/// \ingroup PkgProperty_map
/// Property map that accesses the first item of a `std::pair`. 
/// \tparam Pair Instance of `std::pair`. 
/// \cgalModels `LvaluePropertyMap`
///
/// \sa `CGAL::Second_of_pair_property_map<Pair>`
template <typename Pair>
struct First_of_pair_property_map
  : put_get_helper_pass_key_by_reference<typename Pair::first_type&, 
                                         First_of_pair_property_map<Pair> >
{
  typedef Pair key_type; ///< typedef to 'Pair'
  typedef typename Pair::first_type value_type; ///< typedef to `Pair::first_type`
  typedef value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< boost::lvalue_property_map_tag

  /// Access a property map element.
  ///
  /// @tparam pair a key whose first item is accessed
  reference operator[](key_type& pair) const { return reference(pair.first); }

  /// Access a property map element.
  ///
  /// @tparam pair a key whose first item is accessed read-only
  const reference operator[](const key_type& pair) const { return reference(pair.first); }
};
#endif

/// Free function to create a `First_of_pair_property_map` property map. 
///
/// \relates First_of_pair_property_map 
template <class Iter> // Type convertible to key_type
First_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>
make_first_of_pair_property_map(Iter)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  return First_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>();
}

#ifdef CGAL_USE_OLD_PAIR_PROPERTY_MAPS
/// \ingroup PkgProperty_map
/// 
/// Property map that accesses the second item of a `std::pair`. 
/// 
/// \tparam Pair Instance of `std::pair`. 
/// 
/// \cgalModels `LvaluePropertyMap`
/// 
/// \sa `CGAL::First_of_pair_property_map<Pair>`
template <typename Pair>
struct Second_of_pair_property_map
  : public boost::put_get_helper<typename Pair::second_type&,
                                 Second_of_pair_property_map<Pair> >
{
  typedef Pair* key_type; ///< typedef to 'Pair*'
  typedef typename Pair::second_type value_type; ///< typedef to `Pair::second_type`
  typedef value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`

  /// Access a property map element.
  ///
  /// @tparam Iter Type convertible to `key_type`.
  template <class Iter>
  reference operator[](Iter pair) const { return reference(pair->second); }
};
#else
/// \ingroup PkgProperty_map
/// 
/// Property map that accesses the second item of a `std::pair`. 
/// 
/// \tparam Pair Instance of `std::pair`. 
/// 
/// \cgalModels `LvaluePropertyMap`
/// 
/// \sa `CGAL::First_of_pair_property_map<Pair>`
template <typename Pair>
struct Second_of_pair_property_map
  : put_get_helper_pass_key_by_reference<typename Pair::second_type&, 
                                         Second_of_pair_property_map<Pair> >
{
  typedef Pair key_type; ///< typedef to 'Pair'
  typedef typename Pair::second_type value_type; ///< typedef to `Pair::first_type`
  typedef value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< boost::lvalue_property_map_tag

  /// Access a property map element.
  ///
  /// @tparam pair a key whose second item is accessed
  reference operator[](key_type& pair) const { return reference(pair.second); }

  /// Access a property map element.
  ///
  /// @tparam pair a key whose second item is accessed read-only
  const reference operator[](const key_type& pair) const { return reference(pair.second); }
};
#endif
/// Free function to create a Second_of_pair_property_map property map.
///
/// \relates Second_of_pair_property_map 
template <class Iter> // Type convertible to key_type
Second_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>
make_second_of_pair_property_map(Iter)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  return Second_of_pair_property_map<typename CGAL::value_type_traits<Iter>::type>();
}


//=========================================================================

#ifdef CGAL_USE_OLD_PAIR_PROPERTY_MAPS
/// \ingroup PkgProperty_map
/// 
/// Property map that accesses the Nth item of a `boost::tuple`. 
/// 
/// \tparam N Index of the item to access.
/// \tparam Tuple Instance of `boost::tuple`.
/// 
/// \cgalModels `LvaluePropertyMap`
template <int N, typename Tuple>
struct Nth_of_tuple_property_map
  : public boost::put_get_helper<typename boost::tuples::element<N,Tuple>::type&,
                                 Nth_of_tuple_property_map<N,Tuple> >
{
  typedef Tuple* key_type; ///< typedef to 'Tuple*'
  typedef typename boost::tuples::element<N,Tuple>::type value_type; ///< typedef to `boost::tuples::element<N,Tuple>::%type`
  typedef value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`

  /// Access a property map element.
  ///
  /// @tparam Iter Type convertible to `key_type`.
  template <class Iter>
  reference operator[](Iter tuple) const { return (reference) tuple->template get<N>(); }
};
#else
/// \ingroup PkgProperty_map
/// 
/// Property map that accesses the Nth item of a `boost::tuple`. 
/// 
/// \tparam N Index of the item to access.
/// \tparam Tuple Instance of `boost::tuple`.
/// 
/// \cgalModels `LvaluePropertyMap`
template <int N, typename Tuple>
struct Nth_of_tuple_property_map
  : put_get_helper_pass_key_by_reference<typename boost::tuples::element<N,Tuple>::type&,
                                         Nth_of_tuple_property_map<N,Tuple> >
{
  typedef Tuple key_type; ///< typedef to 'Tuple'
  typedef typename boost::tuples::element<N,Tuple>::type value_type; ///< typedef to `boost::tuples::element<N,Tuple>::%type`
  typedef value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`

  /// Access a property map element.
  ///
  /// @tparam tuple a key whose Nth item is accessed
  reference operator[](key_type& tuple) const { return reference(tuple.template get<N>()); }

  /// Access a property map element.
  ///
  /// @tparam tuple a key whose Nth item is accessed read-only
  const reference operator[](const key_type& tuple) const { return reference(tuple.template get<N>()); }
};
#endif

/// Free function to create a Nth_of_tuple_property_map property map.
///
/// \relates Nth_of_tuple_property_map
template <int N, class Iter> // Type convertible to key_type
Nth_of_tuple_property_map<N, typename CGAL::value_type_traits<Iter>::type>
make_nth_of_tuple_property_map(Iter)
{
  // value_type_traits is a workaround as back_insert_iterator's `value_type` is void
  return Nth_of_tuple_property_map<N, typename CGAL::value_type_traits<Iter>::type>();
}


} // namespace CGAL

#endif // CGAL_POINT_SET_PROPERTY_MAP_H

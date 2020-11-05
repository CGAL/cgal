// Copyright (c) 2008-2009 GeometryFactory and INRIA
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
  #include <boost/vector_property_map.hpp>

#endif
#include <boost/tuple/tuple.hpp>
#include <CGAL/tuple.h>

#include <utility> // defines std::pair

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/Iterator_range.h>
#include <CGAL/Cartesian_converter_fwd.h>
#include <CGAL/Kernel_traits_fwd.h>
#include <CGAL/assertions.h>

namespace CGAL {

/// \cond SKIP_DOXYGEN

/// A boolean property map return a const value at compile time
template <typename Key, bool default_value>
class Static_boolean_property_map
{
public:
  typedef Key key_type;
  typedef bool value_type;
  typedef bool reference;
  typedef boost::read_write_property_map_tag category;

public:

  inline friend
  value_type
  get(Static_boolean_property_map, const key_type&)
  {
    return default_value;
  }

  inline friend
  void
  put(Static_boolean_property_map, const key_type&, const value_type&)
  {}
};


template <typename PM1, typename PM2>
class OR_property_map {
  typedef typename PM1::key_type key_type;
  typedef typename PM1::value_type value_type;
  typedef typename PM1::reference reference;
  typedef boost::read_write_property_map_tag category;

  PM1 pm1;
  PM2 pm2;

 public:
  OR_property_map() {} // required by boost::connected_components

  OR_property_map(PM1 pm1, PM2 pm2)
    : pm1(pm1),pm2(pm2)
  {}

  inline friend
    value_type
    get(const OR_property_map& pm, const key_type& k)
  {
    return get(pm.pm1,k) || get(pm.pm2,k);
  }

  inline friend
    void
  put(OR_property_map& pm, const key_type& k, const value_type& v)
  {
    put(pm.pm1,k, v);
    put(pm.pm2,k, v);
  }
};

template <class PM1, class PM2>
OR_property_map<PM1, PM2>
make_OR_property_map(const PM1& pm1, const PM2& pm2)
{
  return OR_property_map<PM1, PM2>(pm1, pm2);
}

// A property map that uses the result of a property map as key.
template <class KeyMap, class ValueMap>
struct Property_map_binder{
  typedef typename boost::property_traits<KeyMap>::key_type key_type;
  typedef typename boost::property_traits<ValueMap>::value_type value_type;
  typedef typename boost::property_traits<ValueMap>::reference reference;
  typedef boost::read_write_property_map_tag category;

  KeyMap key_map;
  ValueMap value_map;

  Property_map_binder(const KeyMap& key_map, const ValueMap& value_map)
    : key_map(key_map)
    , value_map(value_map)
  {}

  friend
  reference get(const Property_map_binder& map, key_type k)
  {
    return get(map.value_map, get(map.key_map,k));
  }
  friend
  void put(const Property_map_binder& map, key_type k, const value_type& v)
  {
    put(map.value_map, get(map.key_map,k), v);
  }
};

template <class KeyMap, class ValueMap>
Property_map_binder<KeyMap, ValueMap>
bind_property_maps(const KeyMap& src, const ValueMap& tgt)
{
  return Property_map_binder<KeyMap, ValueMap>(src, tgt);
}


/// Property map that accesses a value from an iterator
///
/// \cgalModels `ReadablePropertyMap`
///
/// \tparam InputIterator an input iterator
/// \endcond
template<class InputIterator>
struct Input_iterator_property_map{
  typedef InputIterator key_type;
  typedef typename std::iterator_traits<InputIterator>::value_type value_type;
  typedef typename std::iterator_traits<InputIterator>::reference reference;
  typedef boost::readable_property_map_tag category;

  /// Free function to use a get the value from an iterator using Input_iterator_property_map.
  inline friend
  reference
  get(Input_iterator_property_map<InputIterator>,InputIterator it){ return *it; }
};

/// \ingroup PkgPropertyMapRef
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
  typedef const value_type& reference; ///< typedef to 'T&'
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`

  /// Access a property map element.
  ///
  /// @tparam Iter Type convertible to `key_type`.
  template <class Iter>
  value_type& operator[](Iter it) const { return reference(*it); }
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

/// \ingroup PkgPropertyMapRef
/// A `LvaluePropertyMap` property map mapping a key to itself (by reference).
///
/// \cgalModels `LvaluePropertyMap`
template <typename T>
struct Identity_property_map
{
  typedef T key_type; ///< typedef to `T`
  typedef T value_type; ///< typedef to `T`
  typedef const T& reference; ///< typedef to `T&`
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`
  /// Access a property map element.
  /// @param k a key which is returned as mapped value.
  value_type& operator[](key_type& k) const { return k; }

  typedef Identity_property_map<T> Self;
  /// \name Put/get free functions
  /// @{
  friend reference get(const Self&,const key_type& k) {return k;}
  friend void put(const Self&,key_type& k, const value_type& v) {k=v;}
  /// @}
};


/// \cond SKIP_IN_MANUAL
template <typename T>
struct Identity_property_map_no_lvalue
{
  typedef T key_type; ///< typedef to `T`
  typedef T value_type; ///< typedef to `T`
  typedef T reference; ///< typedef to `T`
  typedef boost::readable_property_map_tag category; ///< `boost::readable_property_map_tag`

  typedef Identity_property_map_no_lvalue<T> Self;

  friend reference get(const Self&, const key_type& k) {return k;}
};
/// \endcond

/// Free function to create a `Identity_property_map` property map.
///
/// \relates Identity_property_map
template <class T> // Key and value type
Identity_property_map<T>
  make_identity_property_map(T)
{
  return Identity_property_map<T>();
}


/// \ingroup PkgPropertyMapRef
/// Property map that accesses the first item of a `std::pair`.
/// \tparam Pair Instance of `std::pair`.
/// \cgalModels `LvaluePropertyMap`
///
/// \sa `CGAL::Second_of_pair_property_map<Pair>`
template <typename Pair>
struct First_of_pair_property_map
{
  typedef Pair key_type; ///< typedef to `Pair`
  typedef typename Pair::first_type value_type; ///< typedef to `Pair::first_type`
  typedef const value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< boost::lvalue_property_map_tag

  /// Access a property map element.
  /// @param pair a key whose first item is accessed
  value_type& operator[](key_type& pair) const { return pair.first; }

  typedef First_of_pair_property_map<Pair> Self;
  /// \name Put/get free functions
  /// @{
  friend reference get(const Self&,const key_type& k) {return k.first;}
  friend void put(const Self&,key_type& k, const value_type& v) {k.first=v;}
  /// @}
};

/// Free function to create a `First_of_pair_property_map` property map.
///
/// \relates First_of_pair_property_map
template <class Pair> // Pair type
First_of_pair_property_map<Pair>
  make_first_of_pair_property_map(Pair)
{
  return First_of_pair_property_map<Pair>();
}

/// \ingroup PkgPropertyMapRef
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
{
  typedef Pair key_type; ///< typedef to `Pair`
  typedef typename Pair::second_type value_type; ///< typedef to `Pair::second_type`
  typedef const value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< boost::lvalue_property_map_tag

  /// Access a property map element.
  /// @param pair a key whose second item is accessed
  value_type& operator[](key_type& pair) const { return pair.second; }

  typedef Second_of_pair_property_map<Pair> Self;
  /// \name Put/get free functions
  /// @{
  friend reference get(const Self&,const key_type& k) {return k.second;}
  friend void put(const Self&,key_type& k, const value_type& v) {k.second=v;}
  /// @}
};

/// Free function to create a Second_of_pair_property_map property map.
///
/// \relates Second_of_pair_property_map
template <class Pair> // Pair type
Second_of_pair_property_map<Pair>
  make_second_of_pair_property_map(Pair)
{
  return Second_of_pair_property_map<Pair>();
}

/// \ingroup PkgPropertyMapRef
///
/// Property map that accesses the Nth item of a `boost::tuple` or a `std::tuple`.
///
/// \tparam N %Index of the item to access.
/// \tparam Tuple Instance of `boost::tuple` or `std::tuple`.
///
/// \cgalModels `LvaluePropertyMap`
template <int N, typename Tuple>
struct Nth_of_tuple_property_map
{
  typedef Tuple key_type; ///< typedef to `Tuple`
  #ifdef DOXYGEN_RUNNING
  typedef unspecified_type value_type;  ///< typedef to the N-th type of the tuple
  #else
  typedef typename boost::tuples::element<N,Tuple>::type value_type;
  #endif
  typedef const value_type& reference; ///< typedef to `value_type&`
  typedef boost::lvalue_property_map_tag category; ///< `boost::lvalue_property_map_tag`
  /// Access a property map element.
  /// @param tuple a key whose Nth item is accessed
  value_type& operator[](key_type& tuple) const { return tuple.template get<N>(); }

  typedef Nth_of_tuple_property_map<N,Tuple> Self;
  /// \name Put/get free functions
  /// @{
  friend reference get(const Self&,const key_type& k) {return k.template get<N>();}
  friend void put(const Self&,key_type& k, const value_type& v) {k.template get<N>()=v;}
  /// @}
};

template <int N, typename ... T>
struct Nth_of_tuple_property_map<N,std::tuple<T...> >
{
  typedef std::tuple<T...> Tuple;
  typedef Tuple key_type;
  typedef typename std::tuple_element<N,Tuple>::type value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  value_type& operator[](key_type& tuple) const { return get<N>(tuple); }

  typedef Nth_of_tuple_property_map<N,Tuple> Self;
  friend reference get(const Self&,const key_type& k) {return std::get<N>(k);}
  friend void put(const Self&,key_type& k, const value_type& v) {std::get<N>(k)=v;}
};

/// Free function to create a Nth_of_tuple_property_map property map.
///
/// \relates Nth_of_tuple_property_map
template <int N, class Tuple> // Tuple type
Nth_of_tuple_property_map<N, Tuple>
  make_nth_of_tuple_property_map(Tuple)
{
  return Nth_of_tuple_property_map<N, Tuple>();
}

/// \ingroup PkgPropertyMapRef
/// Struct that turns a property map into a unary functor with
/// `operator()(key k)` calling the get function with `k`
template <class PropertyMap>
struct Property_map_to_unary_function{
  typedef typename boost::property_traits<PropertyMap>::key_type argument_type;
  typedef typename boost::property_traits<PropertyMap>::reference result_type;

  PropertyMap map;
  Property_map_to_unary_function(PropertyMap m=PropertyMap())
    : map(m)
  {}

  template <class KeyType>
  result_type
  operator()(const KeyType& a) const
  {
    return get(map,a);
  }
};

/// \ingroup PkgPropertyMapRef
/// Utility class providing shortcuts to property maps based on raw pointers
template <class T>
struct Pointer_property_map{
  typedef boost::iterator_property_map< T*,
                              boost::typed_identity_property_map<std::size_t>,
                              T,
                              T&> type; ///< mutable `LvaluePropertyMap`
  typedef boost::iterator_property_map< const T*,
                              boost::typed_identity_property_map<std::size_t>,
                              T,
                              const T&> const_type; ///< non-mutable `LvaluePropertyMap`
};

/// \ingroup PkgPropertyMapRef
/// Starting from boost 1.55, the use of raw pointers as property maps has been deprecated.
/// This function is a shortcut to the recommanded replacement:
/// `boost::make_iterator_property_map(<pointer>, boost::typed_identity_property_map<std::size_t>())`
/// Note that the property map is a mutable `LvaluePropertyMap` with `std::size_t` as key.
template <class T>
inline
typename Pointer_property_map<T>::type
make_property_map(T* pointer)
{
  return typename Pointer_property_map<T>::type(pointer);
}

/// \ingroup PkgPropertyMapRef
/// equivalent to `make_property_map(&v[0])`
/// Note that `v` must not be modified while using the property map created
template <class T>
inline
typename Pointer_property_map<T>::type
make_property_map(std::vector<T>& v)
{
  if(v.empty()){
    return make_property_map(static_cast<T*>(nullptr));
  }
  return make_property_map(&v[0]);
}

/// \ingroup PkgPropertyMapRef
/// Non-mutable version
template <class T>
inline
typename Pointer_property_map<T>::const_type
make_property_map(const T* pointer)
{
  return typename Pointer_property_map<T>::const_type(pointer);
}

/// \ingroup PkgPropertyMapRef
/// equivalent to `make_property_map(&v[0])`
/// Note that `v` must not be modified while using the property map created
template <class T>
inline
typename Pointer_property_map<T>::const_type
make_property_map(const std::vector<T>& v)
{
  return make_property_map(&v[0]);
}

/// \ingroup PkgPropertyMapRef
/// Property map that returns a fixed value.
/// Note that this value is chosen when the map is constructed and cannot
/// be changed afterwards. Specifically, the free function `put()` does nothing.
///
/// \cgalModels `ReadWritePropertyMap`
template<class KeyType, class ValueType>
struct Constant_property_map
{
  ValueType default_value;

  typedef KeyType                                       key_type;
  typedef ValueType                                     value_type;
  typedef value_type&                                   reference;
  typedef boost::read_write_property_map_tag            category;

  Constant_property_map(const value_type& default_value = value_type()) : default_value (default_value) { }

  /// Free function that returns `pm.default_value`.
  inline friend value_type
  get (const Constant_property_map& pm, const key_type&){ return pm.default_value; }

  /// Free function that does nothing.
  inline friend void
  put (const Constant_property_map&, const key_type&, const value_type&) { }
};

/// \ingroup PkgPropertyMapRef
/// Read-write property map turning a set (such a `std::set`,
/// `boost::unordered_set`, `std::unordered_set`) into a property map
/// associating a Boolean to the value type of the set. The function `get` will
/// return `true` if the key is inside the set and `false` otherwise. The `put`
/// function will insert an element in the set if `true` is passed and erase it
/// otherwise.
/// \cgalModels `ReadWritePropertyMap`
template<class Set>
struct Boolean_property_map
{
  typedef typename Set::value_type key_type;
  typedef bool value_type;
  typedef bool reference;
  typedef boost::read_write_property_map_tag category;

  Set* set_ptr;
  /// Constructor taking a copy of the set. Note that `set_` must be valid
  /// while the property map is in use.
  Boolean_property_map(Set& set_) : set_ptr(&set_) {}
  Boolean_property_map() : set_ptr(nullptr) {}

  friend bool get(const Boolean_property_map<Set>& pm, const key_type& k)
  {
    CGAL_assertion(pm.set_ptr!=nullptr);
    return pm.set_ptr->count(k) != 0;
  }

  friend void put(Boolean_property_map<Set>& pm, const key_type& k, bool v)
  {
    CGAL_assertion(pm.set_ptr!=nullptr);
    if (v)
      pm.set_ptr->insert(k);
    else
      pm.set_ptr->erase(k);
  }
};

/// \ingroup PkgPropertyMapRef
/// returns `Boolean_property_map<Set>(set_)`
template <class Set>
Boolean_property_map<Set>
make_boolean_property_map(Set& set_)
{
  return Boolean_property_map<Set>(set_);
}

/// \ingroup PkgPropertyMapRef
/// Read-write property map doing on-the-fly conversions between two default constructible \cgal %Cartesian kernels.
/// Its value type is `GeomObject` and its key type is the same as `Vpm`.
/// `GeomObject` must be a geometric object from a \cgal kernel.
/// `Vpm` is a model `of ReadWritePropertyMap` and its value type must be
/// a geometric object of the same type as `GeomObject` but possibly from
/// another kernel.
/// Conversions between the two geometric objects are done using `Cartesian_converter`.
/// \cgalModels `ReadWritePropertyMap`
template<class GeomObject, class Vpm>
struct Cartesian_converter_property_map
{
  typedef typename boost::property_traits<Vpm>::key_type key_type;
  typedef GeomObject value_type;
  typedef value_type reference;
  typedef boost::read_write_property_map_tag category;
  Vpm vpm;

  typedef typename Kernel_traits<GeomObject>::type K2;
  typedef typename Kernel_traits<typename boost::property_traits<Vpm>::value_type>::type K1;

  Cartesian_converter_property_map(Vpm vpm):vpm(vpm){}

  friend value_type get(const Cartesian_converter_property_map<GeomObject, Vpm>& pm, const key_type& k)
  {
    return
     CGAL::Cartesian_converter<K1, K2>()(get(pm.vpm, k));
  }

  friend void put(Cartesian_converter_property_map<GeomObject, Vpm>& pm, const key_type& k, const value_type& v)
  {
    put(pm.vpm, k, CGAL::Cartesian_converter<K2, K1>()(v));
  }
};

/// \ingroup PkgPropertyMapRef
/// returns `Cartesian_converter_property_map<GeomObject, Vpm>(vpm)`
template<class GeomObject, class Vpm>
Cartesian_converter_property_map<GeomObject, Vpm>
make_cartesian_converter_property_map(Vpm vpm)
{
  return Cartesian_converter_property_map<GeomObject, Vpm>(vpm);
}

/// \cond SKIP_IN_MANUAL
// Syntaxic sugar for transform_iterator+pmap_to_unary_function
template <typename Iterator, typename Pmap>
typename boost::transform_iterator<CGAL::Property_map_to_unary_function<Pmap>, Iterator>
make_transform_iterator_from_property_map (Iterator it, Pmap pmap)
{
  return boost::make_transform_iterator (it, CGAL::Property_map_to_unary_function<Pmap>(pmap));
}

// Syntaxic sugar for make_range+transform_iterator+pmap_to_unary_function
template <typename Range, typename Pmap>
CGAL::Iterator_range<typename boost::transform_iterator<CGAL::Property_map_to_unary_function<Pmap>,
                                                        typename Range::const_iterator> >
make_transform_range_from_property_map (const Range& range, Pmap pmap)
{
  return CGAL::make_range
    (make_transform_iterator_from_property_map (range.begin(), pmap),
     make_transform_iterator_from_property_map (range.end(), pmap));
}

// Syntaxic sugar for make_range+transform_iterator+pmap_to_unary_function
template <typename Range, typename Pmap>
CGAL::Iterator_range<typename boost::transform_iterator<CGAL::Property_map_to_unary_function<Pmap>,
                                                        typename Range::iterator> >
make_transform_range_from_property_map (Range& range, Pmap pmap)
{
  return CGAL::make_range
    (make_transform_iterator_from_property_map (range.begin(), pmap),
     make_transform_iterator_from_property_map (range.end(), pmap));
}

template <typename SizeType>
CGAL::Iterator_range<boost::counting_iterator<SizeType> >
make_counting_range (const SizeType begin, const SizeType end)
{
  return CGAL::make_range (boost::counting_iterator<SizeType>(begin),
                           boost::counting_iterator<SizeType>(end));
}

/// \endcond

} // namespace CGAL



#endif // CGAL_POINT_SET_PROPERTY_MAP_H

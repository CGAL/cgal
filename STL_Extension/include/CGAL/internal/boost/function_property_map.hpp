//
//=======================================================================
// Author: Philipp Moeller
//
// Copyright 2012, Philipp Moeller
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//

// Note: This file was copied from Boost v1.64 because it was introduced
//       in a release (~v1.51) that is newer than the oldest release of Boost
//       that CGAL supports.
//       the property map 'function_property_map' is (currently) used
//       in the packages Triangulation_2 and Triangulation_3.

#ifndef CGAL_INTERNAL_BOOST_PROPERTY_MAP_FUNCTION_PROPERTY_MAP_HPP
#define CGAL_INTERNAL_BOOST_PROPERTY_MAP_FUNCTION_PROPERTY_MAP_HPP

#include <boost/config.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <utility>

namespace CGAL {
namespace internal {
namespace boost_ {

template<typename Func,
         typename Key,
         typename Ret = typename boost::result_of<const Func(const Key&)>::type>
class function_property_map
  : public boost::put_get_helper<Ret, function_property_map<Func, Key, Ret> >
{
public:
  typedef Key key_type;
  typedef Ret reference;
  typedef typename boost::remove_cv<typename boost::remove_reference<Ret>::type>::type value_type;

  typedef typename boost::mpl::if_<
                     boost::mpl::and_<
                       boost::is_reference<Ret>,
                       boost::mpl::not_< boost::is_const<Ret> >
                     >,
                     boost::lvalue_property_map_tag,
                     boost::readable_property_map_tag>::type
    category;

  function_property_map(Func f = Func()) : f(f) {}

  reference operator[](const Key& k) const {
    return f(k);
  }

private:
  Func f;
};

template<typename Key, typename Func>
function_property_map<Func, Key>
make_function_property_map(const Func& f) {
  return function_property_map<Func, Key>(f);
}

template<typename Key, typename Ret, typename Func>
function_property_map<Func, Key, Ret>
make_function_property_map(const Func& f) {
  return function_property_map<Func, Key, Ret>(f);
}

} // boost_
} // internal
} // CGAL

#endif /* CGAL_INTERNAL_BOOST_PROPERTY_MAP_FUNCTION_PROPERTY_MAP_HPP */

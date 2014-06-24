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

#ifndef CGAL_PROPERTY_MAP_FUNCTION_PROPERTY_MAP_H
#define CGAL_PROPERTY_MAP_FUNCTION_PROPERTY_MAP_H

// This is only in boost from 1.50 onward, so it is included here.  It
// has a different include guard and is in a different namespace as
// well to prevent it from affecting users.

#include <boost/config.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>

#include <CGAL/property_map.h>

#include <utility>

namespace CGAL {

template<typename Func, typename Key, typename Ret = typename boost::result_of<const Func(const Key&)>::type>
class function_property_map: public boost::put_get_helper<Ret, function_property_map<Func, Key, Ret> > {
  public:
  typedef Key key_type;
  typedef Ret reference;
  typedef typename boost::remove_cv<typename boost::remove_reference<Ret>::type>::type value_type;

  typedef typename boost::mpl::if_<
                     boost::mpl::and_<
                       boost::is_reference<Ret>,
                       boost::mpl::not_<boost::is_const<Ret> >
                     >,
                     boost::lvalue_property_map_tag,
                     boost::readable_property_map_tag>::type
    category;

  function_property_map(Func f = Func()) : f(f) {}

  reference operator[](key_type k) const {
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

} // CGAL

#endif /* CGAL_PROPERTY_MAP_FUNCTION_PROPERTY_MAP_H */

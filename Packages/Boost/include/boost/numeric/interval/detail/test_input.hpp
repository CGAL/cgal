/* Boost interval/detail/test_input.hpp file
 *
 * Copyright Hervé Brönnimann, Guillaume Melquiond, Sylvain Pion 2002
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation,
 *
 * None of the above authors nor Polytechnic University make any
 * representation about the suitability of this software for any
 * purpose. It is provided "as is" without express or implied warranty.
 *
 * $Id$
 */

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_TEST_INPUT_HPP
#define BOOST_NUMERIC_INTERVAL_DETAIL_TEST_INPUT_HPP

#include <boost/numeric/interval/detail/interval_prototype.hpp>

namespace boost {
namespace numeric {
namespace interval_lib {
namespace detail {

template<class T> inline
bool is_neg(const T& x) { return x < static_cast<T>(0); }

template<class T> inline
bool is_zero(const T& x) { return x == static_cast<T>(0); }

template<class T> inline
bool is_pos(const T& x) { return x > static_cast<T>(0); }

template<class T, class Policies> inline
bool test_input(const interval<T, Policies>& x) {
  typedef typename Policies::checking checking;
  return checking::is_empty(x.lower(), x.upper());
}

template<class T, class Policies1, class Policies2> inline
bool test_input(const interval<T, Policies1>& x, const interval<T, Policies2>& y) {
  typedef typename Policies1::checking checking1;
  typedef typename Policies2::checking checking2;
  return checking1::is_empty(x.lower(), x.upper()) ||
         checking2::is_empty(y.lower(), y.upper());
}

template<class T, class Policies> inline
bool test_input(const T& x, const interval<T, Policies>& y) {
  typedef typename Policies::checking checking;
  return checking::is_nan(x) || checking::is_empty(y.lower(), y.upper());
}

template<class T, class Policies> inline
bool test_input(const interval<T, Policies>& x, const T& y) {
  typedef typename Policies::checking checking;
  return checking::is_empty(x.lower(), x.upper()) || checking::is_nan(y);
}

template<class T, class Policies> inline
bool test_input(const T& x) {
  typedef typename Policies::checking checking;
  return checking::is_nan(x);
}

template<class T, class Policies> inline
bool test_input(const T& x, const T& y) {
  typedef typename Policies::checking checking;
  return checking::is_nan(x) || checking::is_nan(y);
}

} // namespace detail
} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_DETAIL_TEST_INPUT_HPP

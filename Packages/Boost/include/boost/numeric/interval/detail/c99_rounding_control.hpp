/* Boost interval/detail/c99_rounding_control.hpp file
 *
 * Copyright Jens Maurer 2000
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

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP
#define BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP

#include <boost/numeric/interval/detail/c99sub_rounding_control.hpp>

namespace boost {
namespace numeric {
namespace interval_lib {
namespace detail {

struct c99_rounding_control
{
  template<class T>
  static T force_rounding(const T& r) { volatile T r_ = r; return r_; }
};

} // namespace detail

template<>
struct rounding_control<float>:
  detail::c99_rounding_control { };

template<>
struct rounding_control<double>:
  detail::c99_rounding_control { };

template<>
struct rounding_control<long double>:
  detail::c99_rounding_control { };

} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP

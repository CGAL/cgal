/* Boost interval/detail/c99sub_rounding_control.hpp file
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

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_C99SUB_ROUNDING_CONTROL_HPP
#define BOOST_NUMERIC_INTERVAL_DETAIL_C99SUB_ROUNDING_CONTROL_HPP

#include <fenv.h>      // ISO C 99 rounding mode control

namespace boost {
namespace numeric {
namespace interval_lib {
namespace detail {

extern "C" { double rint(double); }

struct c99_rounding
{
  typedef int rounding_mode;

  static void set_rounding_mode(const rounding_mode mode) { fesetround(mode); }
  static void get_rounding_mode(rounding_mode &mode) { mode = fegetround(); }
  static void downward()    { set_rounding_mode(FE_DOWNWARD);   }
  static void upward()      { set_rounding_mode(FE_UPWARD);     }
  static void to_nearest()  { set_rounding_mode(FE_TONEAREST);  }
  static void toward_zero() { set_rounding_mode(FE_TOWARDZERO); }

  template<class T>
  static T to_int(const T& r) { return rint(r); }
};

} // namespace detail
} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_DETAIL_C99SUB_ROUBDING_CONTROL_HPP

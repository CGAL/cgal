/* Boost interval/detail/x86gcc_rounding_control.hpp file
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

#ifndef BOOST_NUMERIC_INTERVAL_EXT_X86_FAST_ROUNDING_CONTROL_HPP
#define BOOST_NUMERIC_INTERVAL_EXT_X86_FAST_ROUNDING_CONTROL_HPP

namespace boost {
namespace numeric {
namespace interval_lib {

namespace detail {

// exceptions masked, expected precision (the mask is 0x0300)
static const fpu_rounding_modes rnd_mode_f = { 0x107f, 0x147f, 0x187f, 0x1c7f };
static const fpu_rounding_modes rnd_mode_d = { 0x127f, 0x167f, 0x1a7f, 0x1e7f };
static const fpu_rounding_modes rnd_mode_l = { 0x137f, 0x177f, 0x1b7f, 0x1f7f };

} // namespace detail

template<class T>
struct x86_fast_rounding_control;

template<>
struct x86_fast_rounding_control<float>: detail::x86_rounding
{
  static void to_nearest()  { set_rounding_mode(detail::rnd_mode_f.to_nearest);  }
  static void downward()    { set_rounding_mode(detail::rnd_mode_f.downward);    }
  static void upward()      { set_rounding_mode(detail::rnd_mode_f.upward);      }
  static void toward_zero() { set_rounding_mode(detail::rnd_mode_f.toward_zero); }
  static const float& force_rounding(const float& r) { return r; }
};

template<>
struct x86_fast_rounding_control<double>: detail::x86_rounding
{
  static void to_nearest()  { set_rounding_mode(detail::rnd_mode_d.to_nearest);  }
  static void downward()    { set_rounding_mode(detail::rnd_mode_d.downward);    }
  static void upward()      { set_rounding_mode(detail::rnd_mode_d.upward);      }
  static void toward_zero() { set_rounding_mode(detail::rnd_mode_d.toward_zero); }
  static const double& force_rounding(const double& r) { return r; }
};

template<>
struct x86_fast_rounding_control<long double>: detail::x86_rounding
{
  static void to_nearest()  { set_rounding_mode(detail::rnd_mode_l.to_nearest);  }
  static void downward()    { set_rounding_mode(detail::rnd_mode_l.downward);    }
  static void upward()      { set_rounding_mode(detail::rnd_mode_l.upward);      }
  static void toward_zero() { set_rounding_mode(detail::rnd_mode_l.toward_zero); }
  static const long double& force_rounding(const long double& r) { return r; }
};

} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_EXT_X86_FAST_ROUNDING_CONTROL_HPP

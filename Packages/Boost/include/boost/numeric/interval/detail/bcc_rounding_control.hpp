/* Boost interval/detail/bcc_rounding_control.hpp file
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

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_BCC_ROUNDING_CONTROL_HPP
#define BOOST_NUMERIC_INTERVAL_DETAIL_BCC_ROUNDING_CONTROL_HPP

#ifndef __BORLANDC__
#  error This header is only intended for Borland C++.
#endif

#ifndef _M_IX86
#  error This header only works on x86 CPUs.
#endif

#include <float.h>      // Borland C++ rounding control

namespace boost {
namespace numeric {
namespace interval_lib {
namespace detail {

#ifndef BOOST_NUMERIC_INTERVAL_KEEP_EXCEPTIONS_FOR_BCC
extern "C" { unsigned int _RTLENTRY _fm_init(void); }

struct borland_workaround {
  borland_workaround() { _fm_init(); }
};

static borland_workaround borland_workaround_exec;
#endif // BOOST_NUMERIC_INTERVAL_KEEP_EXCEPTIONS_FOR_BCC

__inline double rint(double)
{ __emit__(0xD9); __emit__(0xFC); /* asm FRNDINT */ }

struct x86_rounding
{
  typedef unsigned int rounding_mode;
  static void get_rounding_mode(rounding_mode& mode)
  { mode = _control87(0, 0); }
  static void set_rounding_mode(const rounding_mode mode)
  { _control87(mode, 0xffff); } 
  static double to_int(const double& x) { return rint(x); }
};

} // namespace detail
} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif /* BOOST_NUMERIC_INTERVAL_DETAIL_BCC_ROUNDING_CONTROL_HPP */

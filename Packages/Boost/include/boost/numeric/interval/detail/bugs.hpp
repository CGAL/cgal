/* Boost interval/detail/bugs.hpp file
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

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_BUGS
#define BOOST_NUMERIC_INTERVAL_DETAIL_BUGS

#include <boost/config.hpp>

#if defined(__GLIBC__) && !defined(__GLIBCPP__) && (defined(__USE_MISC) || defined(__USE_XOPEN_EXTENDED) || defined(__USE_ISOC99)) && !defined(__ICC)
#  define BOOST_HAS_INV_HYPERBOLIC
#endif

#ifdef BOOST_NO_STDC_NAMESPACE
#  define BOOST_NUMERIC_INTERVAL_using_math(a) using ::a
#  ifdef BOOST_HAS_INV_HYPERBOLIC
#    define BOOST_NUMERIC_INTERVAL_using_ahyp(a) using ::a
#  endif
#else
#  define BOOST_NUMERIC_INTERVAL_using_math(a) using std::a
#  if defined(__GNUC__) && (__GNUC__ == 3) && (__GNUC_MINOR__ == 4)
#    define BOOST_NUMERIC_INTERVAL_using_ahyp(a) using ::a
#  elif defined(BOOST_HAS_INV_HYPERBOLIC)
#    define BOOST_NUMERIC_INTERVAL_using_ahyp(a) using std::a
#  endif
#endif

#if defined(__COMO__) || defined(BOOST_INTEL)
#  define BOOST_NUMERIC_INTERVAL_using_max(a) using std::a
#elif defined(BOOST_NO_STDC_NAMESPACE)
#  define BOOST_NUMERIC_INTERVAL_using_max(a) using ::a
#else
#  define BOOST_NUMERIC_INTERVAL_using_max(a) using std::a
#endif

#ifndef BOOST_NUMERIC_INTERVAL_using_ahyp
#  define BOOST_NUMERIC_INTERVAL_using_ahyp(a)
#endif

#if defined(__GNUC__) && (__GNUC__ <= 2)
// cf PR c++/1981 for a description of the bug
#include <algorithm>
#include <cmath>
namespace boost {
namespace numeric {
  using std::min;
  using std::max;
  using std::sqrt;
  using std::exp;
  using std::log;
  using std::cos;
  using std::tan;
  using std::asin;
  using std::acos;
  using std::atan;
  using std::ceil;
  using std::floor;
  using std::sinh;
  using std::cosh;
  using std::tanh;
# undef BOOST_NUMERIC_INTERVAL_using_max
# undef BOOST_NUMERIC_INTERVAL_using_math
# define BOOST_NUMERIC_INTERVAL_using_max(a)
# define BOOST_NUMERIC_INTERVAL_using_math(a)
# ifdef BOOST_HAS_INV_HYPERBOLIC
  using std::asinh;
  using std::acosh;
  using std::atanh;
# undef BOOST_NUMERIC_INTERVAL_using_ahyp
# define BOOST_NUMERIC_INTERVAL_using_ahyp(a)
# endif
} // namespace numeric
} // namespace boost
#endif

#endif // BOOST_NUMERIC_INTERVAL_DETAIL_BUGS

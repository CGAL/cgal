/* Boost interval/hw_rounding.hpp template implementation file
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

#ifndef BOOST_NUMERIC_INTERVAL_HW_ROUNDING_HPP
#define BOOST_NUMERIC_INTERVAL_HW_ROUNDING_HPP

#include <boost/numeric/interval/rounding.hpp>
#include <boost/numeric/interval/rounded_arith.hpp>

// define appropriate specialization of rounding_control for built-in types
#if defined(__i386__) || defined(_M_IX86) || defined(__BORLANDC__)
#  include <boost/numeric/interval/detail/x86_rounding_control.hpp>
#elif defined(powerpc) || defined(__powerpc__) || defined(__ppc__)
#  include <boost/numeric/interval/detail/ppc_rounding_control.hpp>
#elif defined(sparc) || defined(__sparc__)
#  include <boost/numeric/interval/detail/sparc_rounding_control.hpp>
#elif defined(__USE_ISOC99)
#  include <boost/numeric/interval/detail/c99_rounding_control.hpp>
#else
#  error Boost.Numeric.Interval: Please specify rounding control mechanism.
#endif

namespace boost {
namespace numeric {
namespace interval_lib {

/*
 * Three specializations of rounded_math<T>
 */

template<>
struct rounded_math<float>
  : save_state<rounded_arith_opp<float> >
{};

template<>
struct rounded_math<double>
  : save_state<rounded_arith_opp<double> >
{};

template<>
struct rounded_math<long double>
  : save_state<rounded_arith_opp<long double> >
{};

} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_HW_ROUNDING_HPP

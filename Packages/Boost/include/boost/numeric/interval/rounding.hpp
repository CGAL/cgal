/* Boost interval/rounding.hpp template implementation file
 *
 * Copyright Hervé Brönnimann, Guillaume Melquiond, Sylvain Pion 2002-2003
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

#ifndef BOOST_NUMERIC_INTERVAL_ROUNDING_HPP
#define BOOST_NUMERIC_INTERVAL_ROUNDING_HPP

namespace boost {
namespace numeric {
namespace interval_lib {

/*
 * Default rounding_control class (does nothing)
 */

template<class T>
struct rounding_control
{
  typedef int rounding_mode;
  static void get_rounding_mode(rounding_mode&) {}
  static void set_rounding_mode(rounding_mode)  {}
  static void upward()     {}
  static void downward()   {}
  static void to_nearest() {}
  static const T& to_int(const T& x)         { return x; }
  static const T& force_rounding(const T& x) { return x; }
};

/*
 * A few rounding control classes (exact/std/opp: see documentation)
 *   rounded_arith_* control the rounding of the arithmetic operators
 *   rounded_transc_* control the rounding of the transcendental functions
 */

template<class T, class Rounding = rounding_control<T> >
struct rounded_arith_exact;

template<class T, class Rounding = rounding_control<T> >
struct rounded_arith_std;

template<class T, class Rounding = rounding_control<T> >
struct rounded_arith_opp;

template<class T, class Rounding>
struct rounded_transc_dummy;

template<class T, class Rounding = rounded_arith_exact<T> > 
struct rounded_transc_exact;

template<class T, class Rounding = rounded_arith_std<T> > 
struct rounded_transc_std;

template<class T, class Rounding = rounded_arith_opp<T> > 
struct rounded_transc_opp;

/*
 * State-saving classes: allow to set and reset rounding control
 */

namespace detail {

template<class Rounding>
struct save_state_unprotected: Rounding
{
  typedef save_state_unprotected<Rounding> unprotected_rounding;
};

} // namespace detail

template<class Rounding>
struct save_state: Rounding
{
  typename Rounding::rounding_mode mode;
  save_state() {
    this->get_rounding_mode(mode);
    this->init();
  }
  ~save_state() { this->set_rounding_mode(mode); }
  typedef detail::save_state_unprotected<Rounding> unprotected_rounding;
};
  
template<class Rounding>
struct save_state_nothing: Rounding
{
  typedef save_state_nothing<Rounding> unprotected_rounding;
};
  
template<class T>
struct rounded_math: save_state_nothing<rounded_arith_exact<T> >
{};

} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_ROUNDING_HPP

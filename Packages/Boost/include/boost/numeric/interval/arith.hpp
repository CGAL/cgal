/* Boost interval/arith.hpp template implementation file
 *
 * Copyright Jens Maurer 2000
 * Copyright Hervé Brönnimann, Guillaume Melquiond, Sylvain Pion 2002-2003
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation.
 *
 * None of the above authors nor Polytechnic University make any
 * representation about the suitability of this software for any
 * purpose. It is provided "as is" without express or implied warranty.
 *
 * $Id$
 */

#ifndef BOOST_NUMERIC_INTERVAL_ARITH_HPP
#define BOOST_NUMERIC_INTERVAL_ARITH_HPP

#include <boost/numeric/interval/detail/interval_prototype.hpp>
#include <boost/numeric/interval/detail/bugs.hpp>
#include <boost/numeric/interval/detail/test_input.hpp>
#include <boost/numeric/interval/detail/division.hpp>
#include <algorithm>

namespace boost {
namespace numeric {

/*
 * Basic arithmetic operators
 */

template<class T, class Policies> inline
const interval<T, Policies>& operator+(const interval<T, Policies>& x)
{
  return x;
}

template<class T, class Policies> inline
interval<T, Policies> operator-(const interval<T, Policies>& x)
{
  if (interval_lib::detail::test_input(x))
    return interval<T, Policies>::empty();
  return interval<T, Policies>(-x.upper(), -x.lower(), true);
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator+=(const interval<T, Policies>& r)
{
  if (interval_lib::detail::test_input(*this, r))
    set_empty();
  else {
    typename Policies::rounding rnd;
    set(rnd.add_down(low, r.low), rnd.add_up(up, r.up));
  }
  return *this;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator+=(const T& r)
{
  if (interval_lib::detail::test_input(*this, r))
    set_empty();
  else {
    typename Policies::rounding rnd;
    set(rnd.add_down(low, r), rnd.add_up(up, r));
  }
  return *this;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator-=(const interval<T, Policies>& r)
{
  if (interval_lib::detail::test_input(*this, r))
    set_empty();
  else {
    typename Policies::rounding rnd;
    set(rnd.sub_down(low, r.up), rnd.sub_up(up, r.low));
  }
  return *this;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator-=(const T& r)
{
  if (interval_lib::detail::test_input(*this, r))
    set_empty();
  else {
    typename Policies::rounding rnd;
    set(rnd.sub_down(low, r), rnd.sub_up(up, r));
  }
  return *this;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator*=(const interval<T, Policies>& r)
{
  return *this = *this * r;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator*=(const T& r)
{
  return *this = r * *this;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator/=(const interval<T, Policies>& r)
{
  return *this = *this / r;
}

template<class T, class Policies> inline
interval<T, Policies>& interval<T, Policies>::operator/=(const T& r)
{
  return *this = *this / r;
}

template<class T, class Policies> inline
interval<T, Policies> operator+(const interval<T, Policies>& x,
                                const interval<T, Policies>& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  typename Policies::rounding rnd;
  return interval<T,Policies>(rnd.add_down(x.lower(), y.lower()),
                              rnd.add_up  (x.upper(), y.upper()), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator+(const T& x, const interval<T, Policies>& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  typename Policies::rounding rnd;
  return interval<T,Policies>(rnd.add_down(x, y.lower()),
                              rnd.add_up  (x, y.upper()), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator+(const interval<T, Policies>& x, const T& y)
{ return y + x; }

template<class T, class Policies> inline
interval<T, Policies> operator-(const interval<T, Policies>& x,
                                const interval<T, Policies>& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  typename Policies::rounding rnd;
  return interval<T,Policies>(rnd.sub_down(x.lower(), y.upper()),
                              rnd.sub_up  (x.upper(), y.lower()), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator-(const T& x, const interval<T, Policies>& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  typename Policies::rounding rnd;
  return interval<T,Policies>(rnd.sub_down(x, y.upper()),
                              rnd.sub_up  (x, y.lower()), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator-(const interval<T, Policies>& x, const T& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  typename Policies::rounding rnd;
  return interval<T,Policies>(rnd.sub_down(x.lower(), y),
                              rnd.sub_up  (x.upper(), y), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator*(const interval<T, Policies>& x,
                                const interval<T, Policies>& y)
{
  BOOST_NUMERIC_INTERVAL_using_max(min);
  BOOST_NUMERIC_INTERVAL_using_max(max);
  typedef interval<T, Policies> I;
  if (interval_lib::detail::test_input(x, y))
    return I::empty();
  typename Policies::rounding rnd;
  const T& xl = x.lower();
  const T& xu = x.upper();
  const T& yl = y.lower();
  const T& yu = y.upper();

  if (interval_lib::detail::is_neg(xl))
    if (interval_lib::detail::is_pos(xu))
      if (interval_lib::detail::is_neg(yl))
        if (interval_lib::detail::is_pos(yu)) // M * M
          return I(min(rnd.mul_down(xl, yu), rnd.mul_down(xu, yl)),
                   max(rnd.mul_up  (xl, yl), rnd.mul_up  (xu, yu)), true);
        else                    // M * N
          return I(rnd.mul_down(xu, yl), rnd.mul_up(xl, yl), true);
      else
        if (interval_lib::detail::is_pos(yu)) // M * P
          return I(rnd.mul_down(xl, yu), rnd.mul_up(xu, yu), true);
        else                    // M * Z
          return I(static_cast<T>(0), static_cast<T>(0), true);
    else
      if (interval_lib::detail::is_neg(yl))
        if (interval_lib::detail::is_pos(yu)) // N * M
          return I(rnd.mul_down(xl, yu), rnd.mul_up(xl, yl), true);
        else                    // N * N
          return I(rnd.mul_down(xu, yu), rnd.mul_up(xl, yl), true);
      else
        if (interval_lib::detail::is_pos(yu)) // N * P
          return I(rnd.mul_down(xl, yu), rnd.mul_up(xu, yl), true);
        else                    // N * Z
          return I(static_cast<T>(0), static_cast<T>(0), true);
  else
    if (interval_lib::detail::is_pos(xu))
      if (interval_lib::detail::is_neg(yl))
        if (interval_lib::detail::is_pos(yu)) // P * M
          return I(rnd.mul_down(xu, yl), rnd.mul_up(xu, yu), true);
        else                    // P * N
          return I(rnd.mul_down(xu, yl), rnd.mul_up(xl, yu), true);
      else
        if (interval_lib::detail::is_pos(yu)) // P * P
          return I(rnd.mul_down(xl, yl), rnd.mul_up(xu, yu), true);
        else                    // P * Z
          return I(static_cast<T>(0), static_cast<T>(0), true);
    else                        // Z * ?
      return I(static_cast<T>(0), static_cast<T>(0), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator*(const T& x, const interval<T, Policies>& y)
{ 
  typedef interval<T, Policies> I;
  if (interval_lib::detail::test_input(x, y))
    return I::empty();
  typename Policies::rounding rnd;
  const T& yl = y.lower();
  const T& yu = y.upper();
  // x is supposed not to be infinite
  if (interval_lib::detail::is_neg(x))
    return I(rnd.mul_down(x, yu), rnd.mul_up(x, yl), true);
  else if (interval_lib::detail::is_zero(x))
    return I(static_cast<T>(0), static_cast<T>(0), true);
  else
    return I(rnd.mul_down(x, yl), rnd.mul_up(x, yu), true);
}

template<class T, class Policies> inline
interval<T, Policies> operator*(const interval<T, Policies>& x, const T& y)
{ return y * x; }

template<class T, class Policies> inline
interval<T, Policies> operator/(const interval<T, Policies>& x,
                                const interval<T, Policies>& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  if (in_zero(y))
    if (!interval_lib::detail::is_zero(y.lower()))
      if (!interval_lib::detail::is_zero(y.upper()))
        return interval_lib::detail::div_zero(x);
      else
        return interval_lib::detail::div_negative(x, y.lower());
    else
      if (!interval_lib::detail::is_zero(y.upper()))
        return interval_lib::detail::div_positive(x, y.upper());
      else
        return interval<T, Policies>::empty();
  else
    return interval_lib::detail::div_non_zero(x, y);
}

template<class T, class Policies> inline
interval<T, Policies> operator/(const T& x, const interval<T, Policies>& y)
{
  if (interval_lib::detail::test_input(x, y))
    return interval<T, Policies>::empty();
  if (in_zero(y))
    if (!interval_lib::detail::is_zero(y.lower()))
      if (!interval_lib::detail::is_zero(y.upper()))
        return interval_lib::detail::div_zero<T, Policies>(x);
      else
        return interval_lib::detail::div_negative<T, Policies>(x, y.lower());
    else
      if (!interval_lib::detail::is_zero(y.upper()))
        return interval_lib::detail::div_positive<T, Policies>(x, y.upper());
      else
        return interval<T, Policies>::empty();
  else
    return interval_lib::detail::div_non_zero(x, y);
}

template<class T, class Policies> inline
interval<T, Policies> operator/(const interval<T, Policies>& x, const T& y)
{
  if (interval_lib::detail::test_input(x, y) || interval_lib::detail::is_zero(y))
    return interval<T, Policies>::empty();
  typename Policies::rounding rnd;
  const T& xl = x.lower();
  const T& xu = x.upper();
  if (interval_lib::detail::is_neg(y))
    return interval<T, Policies>(rnd.div_down(xu, y), rnd.div_up(xl, y), true);
  else
    return interval<T, Policies>(rnd.div_down(xl, y), rnd.div_up(xu, y), true);
}

} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_ARITH_HPP

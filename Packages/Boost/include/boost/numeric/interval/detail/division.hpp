/* Boost interval/detail/division.hpp file
 *
 * Copyright Guillaume Melquiond, Sylvain Pion 2003
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation,
 *
 * None of the above authors make any representation about the
 * suitability of this software for any purpose. It is provided "as
 * is" without express or implied warranty.
 *
 * $Id$
 */

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_DIVISION_HPP
#define BOOST_NUMERIC_INTERVAL_DETAIL_DIVISION_HPP

#include <boost/numeric/interval/detail/interval_prototype.hpp>
#include <boost/numeric/interval/detail/bugs.hpp>
#include <boost/numeric/interval/rounded_arith.hpp>
#include <algorithm>

namespace boost {
namespace numeric {
namespace interval_lib {
namespace detail {

template<class T, class Policies> inline
interval<T, Policies> div_non_zero(const interval<T, Policies>& x,
                                   const interval<T, Policies>& y)
{
  // assert(!in_zero(y));
  typename Policies::rounding rnd;
  typedef interval<T, Policies> I;
  const T& xl = x.lower();
  const T& xu = x.upper();
  const T& yl = y.lower();
  const T& yu = y.upper();
  if (is_neg(xu))
    if (is_neg(yu))
      return I(rnd.div_down(xu, yl), rnd.div_up(xl, yu), true);
    else
      return I(rnd.div_down(xl, yl), rnd.div_up(xu, yu), true);
  else if (is_neg(xl))
    if (is_neg(yu))
      return I(rnd.div_down(xu, yu), rnd.div_up(xl, yu), true);
    else
      return I(rnd.div_down(xl, yl), rnd.div_up(xu, yl), true);
  else
    if (is_neg(yu))
      return I(rnd.div_down(xu, yu), rnd.div_up(xl, yl), true);
    else
      return I(rnd.div_down(xl, yu), rnd.div_up(xu, yl), true);
}

template<class T, class Policies> inline
interval<T, Policies> div_non_zero(const T& x, const interval<T, Policies>& y)
{
  // assert(!in_zero(y));
  typename Policies::rounding rnd;
  typedef interval<T, Policies> I;
  const T& yl = y.lower();
  const T& yu = y.upper();
  if (is_neg(x))
    return I(rnd.div_down(x, yl), rnd.div_up(x, yu), true);
  else
    return I(rnd.div_down(x, yu), rnd.div_up(x, yl), true);
}

template<class T, class Policies> inline
interval<T, Policies> div_positive(const interval<T, Policies>& x, const T& yu)
{
  // assert(yu > T(0));
  if (is_zero(x.lower()) && is_zero(x.upper())) return x;
  typename Policies::rounding rnd;
  typedef interval<T, Policies> I;
  const T& xl = x.lower();
  const T& xu = x.upper();
  typedef typename Policies::checking checking;
  const T& inf = checking::inf();
  if (is_neg(xu))
    return I(-inf, rnd.div_up(xu, yu), true);
  else if (is_neg(xl))
    return I(-inf, inf, true);
  else
    return I(rnd.div_down(xl, yu), inf, true);
}

template<class T, class Policies> inline
interval<T, Policies> div_positive(const T& x, const T& yu)
{
  // assert(yu > T(0));
  typedef interval<T, Policies> I;
  if (is_zero(x)) return I(0, 0, true);
  typename Policies::rounding rnd;
  typedef typename Policies::checking checking;
  const T& inf = checking::inf();
  if (is_neg(x))
    return I(-inf, rnd.div_up(x, yu), true);
  else
    return I(rnd.div_down(x, yu), inf, true);
}

template<class T, class Policies> inline
interval<T, Policies> div_negative(const interval<T, Policies>& x, const T& yl)
{
  // assert(yl < T(0));
  if (is_zero(x.lower()) && is_zero(x.upper()))
    return x;
  typename Policies::rounding rnd;
  typedef interval<T, Policies> I;
  const T& xl = x.lower();
  const T& xu = x.upper();
  typedef typename Policies::checking checking;
  const T& inf = checking::inf();
  if (is_neg(xu))
    return I(rnd.div_down(xu, yl), inf, true);
  else if (is_neg(xl))
    return I(-inf, inf, true);
  else
    return I(-inf, rnd.div_up(xl, yl), true);
}

template<class T, class Policies> inline
interval<T, Policies> div_negative(const T& x, const T& yl)
{
  // assert(yl < T(0));
  typedef interval<T, Policies> I;
  if (is_zero(x)) return I(0, 0, true);
  typename Policies::rounding rnd;
  typedef typename Policies::checking checking;
  const T& inf = checking::inf();
  if (is_neg(x))
    return I(rnd.div_down(x, yl), inf, true);
  else
    return I(-inf, rnd.div_up(x, yl), true);
}

template<class T, class Policies> inline
interval<T, Policies> div_zero(const interval<T, Policies>& x)
{
  if (is_zero(x.lower()) && is_zero(x.upper()))
    return x;
  else return interval<T, Policies>::whole();
}

template<class T, class Policies> inline
interval<T, Policies> div_zero(const T& x)
{
  if (is_zero(x)) return interval<T, Policies>(0, 0, true);
  else return interval<T, Policies>::whole();
}

template<class T, class Policies> inline
interval<T, Policies> div_zero_part1(const interval<T, Policies>& x,
                                     const interval<T, Policies>& y, bool& b)
{
  // assert(y.lower() < 0 && y.upper() > 0);
  if (is_zero(x.lower()) && is_zero(x.upper()))
    { b = false; return x; }
  typename Policies::rounding rnd;
  typedef interval<T, Policies> I;
  const T& xl = x.lower();
  const T& xu = x.upper();
  const T& yl = y.lower();
  const T& yu = y.upper();
  typedef typename Policies::checking checking;
  const T& inf = checking::inf();
  if (is_neg(xu))
    { b = true;  return I(-inf, rnd.div_up(xu, yu), true); }
  else if (is_neg(xl))
    { b = false; return I(-inf, inf, true); }
  else
    { b = true;  return I(-inf, rnd.div_up(xl, yl), true); }
}

template<class T, class Policies> inline
interval<T, Policies> div_zero_part2(const interval<T, Policies>& x,
                                     const interval<T, Policies>& y)
{
  // assert(y.lower() < 0 && y.upper() > 0 && (div_zero_part1(x, y, b), b));
  typename Policies::rounding rnd;
  typedef interval<T, Policies> I;
  typedef typename Policies::checking checking;
  const T& inf = checking::inf();
  if (is_neg(x.upper()))
    return I(rnd.div_down(x.upper(), y.lower()), inf, true);
  else
    return I(rnd.div_down(x.lower(), y.upper()), inf, true);
}

} // namespace detail
} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_DETAIL_DIVISION_HPP

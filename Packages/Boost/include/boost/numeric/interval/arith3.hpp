/* Boost interval/arith3.hpp template implementation file
 *
 * Copyright Guillaume Melquiond 2003
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation.
 *
 * None of the above authors make any representation about the suitability
 * of this software for any purpose. It is provided "as is" without express
 * or implied warranty.
 *
 * $Id$
 */

#ifndef BOOST_NUMERIC_INTERVAL_ARITH3_HPP
#define BOOST_NUMERIC_INTERVAL_ARITH3_HPP

#include <boost/numeric/interval/detail/interval_prototype.hpp>
#include <boost/numeric/interval/detail/test_input.hpp>

namespace boost {
namespace numeric {
namespace interval_lib {

template<class I> inline
I add(const typename I::base_type& x, const typename I::base_type& y)
{
  typedef typename I::traits_type Policies;
  if (detail::test_input<typename I::base_type, Policies>(x, y))
    return I::empty();
  typename Policies::rounding rnd;
  return I(rnd.add_down(x, y), rnd.add_up(x, y), true);
}

template<class I> inline
I sub(const typename I::base_type& x, const typename I::base_type& y)
{
  typedef typename I::traits_type Policies;
  if (detail::test_input<typename I::base_type, Policies>(x, y))
    return I::empty();
  typename Policies::rounding rnd;
  return I(rnd.sub_down(x, y), rnd.sub_up(x, y), true);
}

template<class I> inline
I mul(const typename I::base_type& x, const typename I::base_type& y)
{
  typedef typename I::traits_type Policies;
  if (detail::test_input<typename I::base_type, Policies>(x, y))
    return I::empty();
  typename Policies::rounding rnd;
  return I(rnd.mul_down(x, y), rnd.mul_up(x, y), true);
}

template<class I> inline
I div(const typename I::base_type& x, const typename I::base_type& y)
{
  typedef typename I::traits_type Policies;
  if (detail::test_input<typename I::base_type, Policies>(x, y) || detail::is_zero(y))
    return I::empty();
  typename Policies::rounding rnd;
  return I(rnd.div_down(x, y), rnd.div_up(x, y), true);
}

} // namespace interval_lib
} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_ARITH3_HPP

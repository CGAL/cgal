/* Boost interval/ext/integer.hpp template implementation file
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

#ifndef BOOST_NUMERIC_INTERVAL_EXT_INTEGER_HPP
#define BOOST_NUMERIC_INTERVAL_EXT_INTEGER_HPP

#include <boost/numeric/interval/detail/interval_prototype.hpp>
#include <boost/numeric/interval/detail/test_input.hpp>

namespace boost {
namespace numeric {

template<class T, class Policies> inline
interval<T, Policies> operator+ (const interval<T, Policies>& x, int y)
{
  return x + static_cast<T>(y);
}

template<class T, class Policies> inline
interval<T, Policies> operator+ (int x, const interval<T, Policies>& y)
{
  return static_cast<T>(x) + y;
}

template<class T, class Policies> inline
interval<T, Policies> operator- (const interval<T, Policies>& x, int y)
{
  return x - static_cast<T>(y);
}

template<class T, class Policies> inline
interval<T, Policies> operator- (int x, const interval<T, Policies>& y)
{
  return static_cast<T>(x) - y;
}

template<class T, class Policies> inline
interval<T, Policies> operator* (const interval<T, Policies>& x, int y)
{
  return x * static_cast<T>(y);
}

template<class T, class Policies> inline
interval<T, Policies> operator* (int x, const interval<T, Policies>& y)
{
  return static_cast<T>(x) * y;
}

template<class T, class Policies> inline
interval<T, Policies> operator/ (const interval<T, Policies>& x, int y)
{
  return x / static_cast<T>(y);
}

template<class T, class Policies> inline
interval<T, Policies> operator/ (int x, const interval<T, Policies>& y)
{
  return static_cast<T>(x) / y;
}

} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_EXT_INTEGER_HPP

/* Boost interval/compare/tribool.hpp template implementation file
 *
 * Copyright Guillaume Melquiond 2002-2003
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

#ifndef BOOST_NUMERIC_INTERVAL_COMPARE_TRIBOOL_HPP
#define BOOST_NUMERIC_INTERVAL_COMPARE_TRIBOOL_HPP

#include <boost/numeric/interval/detail/interval_prototype.hpp>
#include <boost/numeric/interval/detail/test_input.hpp>
#include <boost/tribool.hpp>

namespace boost {
namespace numeric {
namespace interval_lib {
namespace compare {
namespace tribool {

template<class T, class Policies1, class Policies2> inline
tribool operator<(const interval<T, Policies1>& x, const interval<T, Policies2>& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() < y.lower()) return true;
  if (x.lower() >= y.upper()) return false;
  return indeterminate;
}

template<class T, class Policies> inline
tribool operator<(const interval<T, Policies>& x, const T& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() < y) return true;
  if (x.lower() >= y) return false;
  return indeterminate;
}

template<class T, class Policies1, class Policies2> inline
tribool operator<=(const interval<T, Policies1>& x, const interval<T, Policies2>& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() <= y.lower()) return true;
  if (x.lower() > y.upper()) return false;
  return indeterminate;
}

template<class T, class Policies> inline
tribool operator<=(const interval<T, Policies>& x, const T& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() <= y) return true;
  if (x.lower() > y) return false;
  return indeterminate;
}

template<class T, class Policies1, class Policies2> inline
tribool operator>(const interval<T, Policies1>& x, const interval<T, Policies2>& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.lower() > y.upper()) return true;
  if (x.upper() <= y.lower()) return false;
  return indeterminate;
}

template<class T, class Policies> inline
tribool operator>(const interval<T, Policies>& x, const T& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.lower() > y) return true;
  if (x.upper() <= y) return false;
  return indeterminate;
}

template<class T, class Policies1, class Policies2> inline
tribool operator>=(const interval<T, Policies1>& x, const interval<T, Policies2>& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.lower() >= y.upper()) return true;
  if (x.upper() < y.lower()) return false;
  return indeterminate;
}

template<class T, class Policies> inline
tribool operator>=(const interval<T, Policies>& x, const T& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.lower() >= y) return true;
  if (x.upper() < y) return false;
  return indeterminate;
}

template<class T, class Policies1, class Policies2> inline
tribool operator==(const interval<T, Policies1>& x, const interval<T, Policies2>& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() == y.lower() && x.lower() == y.upper()) return true;
  if (x.upper() < y.lower() || x.lower() > y.upper()) return false;
  return indeterminate;
}

template<class T, class Policies> inline
tribool operator==(const interval<T, Policies>& x, const T& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() == y && x.lower() == y) return true;
  if (x.upper() < y || x.lower() > y) return false;
  return indeterminate;
}

template<class T, class Policies1, class Policies2> inline
tribool operator!=(const interval<T, Policies1>& x, const interval<T, Policies2>& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() < y.lower() || x.lower() > y.upper()) return true;
  if (x.upper() == y.lower() && x.lower() == y.upper()) return false;
  return indeterminate;
}

template<class T, class Policies> inline
tribool operator!=(const interval<T, Policies>& x, const T& y)
{
  if (detail::test_input(x, y)) throw comparison_error();
  if (x.upper() < y || x.lower() > y) return true;
  if (x.upper() == y && x.lower() == y) return false;
  return indeterminate;
}

} // namespace tribool
} // namespace compare
} // namespace interval_lib
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_INTERVAL_COMPARE_TRIBOOL_HPP

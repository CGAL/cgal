/* boost random/uniform_real.hpp header file
 *
 * Copyright Jens Maurer 2000-2001
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation,
 *
 * Jens Maurer makes no representations about the suitability of this
 * software for any purpose. It is provided "as is" without express or
 * implied warranty.
 *
 * See http://www.boost.org for most recent version including documentation.
 *
 * $Id$
 *
 * Revision history
 *  2001-04-08  added min<max assertion (N. Becker)
 *  2001-02-18  moved to individual header files
 */

#ifndef BOOST_RANDOM_UNIFORM_REAL_HPP
#define BOOST_RANDOM_UNIFORM_REAL_HPP

#include <cassert>
#include <iostream>
#include <boost/config.hpp>
#include <boost/limits.hpp>
#include <boost/static_assert.hpp>

namespace boost {

// uniform distribution on a real range
template<class RealType = double>
class uniform_real
{
public:
  typedef RealType input_type;
  typedef RealType result_type;

  explicit uniform_real(RealType min = RealType(0),
                        RealType max = RealType(1)) 
    : _min(min), _max(max)
  {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif
    assert(min < max);
  }

  // compiler-generated copy ctor and assignment operator are fine

  result_type min() const { return _min; }
  result_type max() const { return _max; }
  void reset() { }

  template<class Engine>
  result_type operator()(Engine& eng) { return eng() * (_max - _min) + _min; }

#if !defined(BOOST_NO_OPERATORS_IN_NAMESPACE) && !defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS)
  template<class CharT, class Traits>
  friend std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os, const uniform_real& ud)
  {
    os << ud._min << " " << ud._max;
    return os;
  }

  template<class CharT, class Traits>
  friend std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT,Traits>& is, uniform_real& ud)
  {
    is >> std::ws >> ud._min >> std::ws >> ud._max;
    return is;
  }
#endif

private:
  RealType _min, _max;
};

} // namespace boost

#endif // BOOST_RANDOM_UNIFORM_REAL_HPP

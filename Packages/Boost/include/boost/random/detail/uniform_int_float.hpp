/* boost random/detail/uniform_int_float.hpp header file
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
 */

#ifndef BOOST_RANDOM_DETAIL_UNIFORM_INT_FLOAT_HPP
#define BOOST_RANDOM_DETAIL_UNIFORM_INT_FLOAT_HPP

#include <boost/random/uniform_01.hpp>


namespace boost {
namespace random {
namespace detail {

template<class UniformRandomNumberGenerator, class IntType = unsigned long>
class uniform_int_float
{
public:
  typedef UniformRandomNumberGenerator base_type;
  typedef IntType result_type;

  uniform_int_float(base_type rng, IntType min = 0, IntType max = 0xffffffff)
    : _rng(rng), _min(min), _max(max)
  {
    init();
  }

  result_type min() const { return _min; }
  result_type max() const { return _max; }
  base_type& base() { return _rng.base(); }
  const base_type& base() const { return _rng.base(); }

  result_type operator()()
  {
    return static_cast<IntType>(_rng() * _range) + _min;
  }

#if !defined(BOOST_NO_OPERATORS_IN_NAMESPACE) && !defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS)
  template<class CharT, class Traits>
  friend std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os, const uniform_int_float& ud)
  {
    os << ud._min << " " << ud._max;
    return os;
  }

  template<class CharT, class Traits>
  friend std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT,Traits>& is, uniform_int_float& ud)
  {
    is >> std::ws >> ud._min >> std::ws >> ud._max;
    ud.init();
    return is;
  }
#endif

private:
  void init()
  {
    _range = static_cast<base_result>(_max-_min)+1;
  }

  typedef typename base_type::result_type base_result;
  uniform_01<base_type> _rng;
  result_type _min, _max;
  base_result _range;
};


} // namespace detail
} // namespace random
} // namespace boost

#endif // BOOST_RANDOM_DETAIL_UNIFORM_INT_FLOAT_HPP

/* Boost interval/limits.hpp template implementation file
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

#ifndef BOOST_NUMERIC_INTERVAL_LIMITS_HPP
#define BOOST_NUMERIC_INTERVAL_LIMITS_HPP

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#include <boost/config.hpp>
#include <boost/limits.hpp>
#include <boost/numeric/interval/detail/interval_prototype.hpp>

namespace std {

template<class T, class Policies>
class numeric_limits<boost::numeric::interval<T, Policies> >
  : public numeric_limits<T>
{
private:
  typedef boost::numeric::interval<T, Policies> I;
  typedef numeric_limits<T> bl;
public:
  static I min() throw() { return I(bl::min(), bl::min()); }
  static I max() throw() { return I(bl::max(), bl::max()); }
  static I epsilon() throw() { return I(bl::epsilon(), bl::epsilon()); }

  BOOST_STATIC_CONSTANT(float_round_style, round_style = round_indeterminate);
  BOOST_STATIC_CONSTANT(bool, is_iec559 = false);

  static I infinity () throw() { return I::whole(); }
  static I quiet_NaN() throw() { return I::empty(); }
  static I signaling_NaN() throw()
  { return I(bl::signaling_NaN(), bl::signaling_Nan()); }
  static I denorm_min() throw()
  { return I(bl::denorm_min(), bl::denorm_min()); }
private:
  static I round_error();    // hide this on purpose, not yet implemented
};

} // namespace std

#endif

#endif // BOOST_NUMERIC_INTERVAL_LIMITS_HPP

/* Boost interval/detail/interval_prototype.hpp file
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

#ifndef BOOST_NUMERIC_INTERVAL_DETAIL_INTERVAL_PROTOTYPE_HPP
#define BOOST_NUMERIC_INTERVAL_DETAIL_INTERVAL_PROTOTYPE_HPP

namespace boost {
namespace numeric {

namespace interval_lib {

template<class T> struct rounded_math;
template<class T> struct checking_strict;
class comparison_error;
template<class Rounding, class Checking> struct policies;

/*
 * default policies class
 */

template<class T>
struct default_policies
{
  typedef policies<rounded_math<T>, checking_strict<T> > type;
};
    
} // namespace interval_lib

template<class T, class Policies = typename interval_lib::default_policies<T>::type >
class interval;  

} // namespace numeric
} // namespace boost

#endif // BOOST_NUMERIC_INTERVAL_DETAIL_INTERVAL_PROTOTYPE_HPP

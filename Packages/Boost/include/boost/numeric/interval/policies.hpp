/* Boost interval/policies.hpp template implementation file
 *
 * Copyright Guillaume Melquiond 2003
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

#ifndef BOOST_NUMERIC_INTERVAL_POLICIES_HPP
#define BOOST_NUMERIC_INTERVAL_POLICIES_HPP

#include <boost/numeric/interval/interval.hpp>

namespace boost {
namespace numeric {
namespace interval_lib {

/*
 * policies class
 */

template<class Rounding, class Checking>
struct policies
{
  typedef Rounding rounding;
  typedef Checking checking;
};

/*
 * policies switching classes
 */

template<class OldInterval, class NewRounding>
class change_rounding
{
  typedef typename OldInterval::base_type T;
  typedef typename OldInterval::traits_type p;
  typedef typename p::checking checking;
public:
  typedef interval<T, policies<NewRounding, checking> > type;
};

template<class OldInterval, class NewChecking>
class change_checking
{
  typedef typename OldInterval::base_type T;
  typedef typename OldInterval::traits_type p;
  typedef typename p::rounding rounding;
public:
  typedef interval<T, policies<rounding, NewChecking> > type;
};

/*
 * Protect / unprotect: control whether the rounding mode is set/reset
 * at each operation, rather than once and for all.
 */

template<class OldInterval>
class unprotect
{
  typedef typename OldInterval::base_type T;
  typedef typename OldInterval::traits_type p;
  typedef typename p::rounding r;
  typedef typename r::unprotected_rounding newRounding;
public:
  typedef typename change_rounding<OldInterval, newRounding>::type type;
};

} // namespace interval_lib
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_INTERVAL_POLICIES_HPP

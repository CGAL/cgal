//-----------------------------------------------------------------------------
// boost variant/detail/substitute_fwd.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_VARIANT_DETAIL_SUBSTITUTE_FWD_HPP
#define BOOST_VARIANT_DETAIL_SUBSTITUTE_FWD_HPP

#include "boost/mpl/aux_/lambda_arity_param.hpp"
#include "boost/mpl/aux_/template_arity.hpp"
#include "boost/mpl/int_fwd.hpp"


///////////////////////////////////////////////////////////////////////////////
// BOOST_VARIANT_DETAIL_NO_SUBSTITUTE
//
// Defined if 'substitute' is not implementable on the current compiler.
//

#include "boost/mpl/aux_/config/ctps.hpp"
#include "boost/mpl/aux_/config/ttp.hpp"

#if defined(BOOST_NO_TEMPLATE_TEMPLATE_PARAMETERS) \
 || defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
 && !defined(BOOST_VARIANT_DETAIL_NO_SUBSTITUTE)
#   define BOOST_VARIANT_DETAIL_NO_SUBSTITUTE
#endif

namespace boost {
namespace detail { namespace variant {

#if !defined(BOOST_VARIANT_DETAIL_NO_SUBSTITUTE)

///////////////////////////////////////////////////////////////////////////////
// metafunction substitute
//
// Substitutes one type for another in the given type expression.
//
template <
      typename T, typename Dest, typename Source
      BOOST_MPL_AUX_LAMBDA_ARITY_PARAM(
          typename Arity = mpl::int_< mpl::aux::template_arity<T>::value >
        )
    >
struct substitute;

#endif // !defined(BOOST_VARIANT_DETAIL_NO_SUBSTITUTE)

}} // namespace detail::variant
} // namespace boost

#endif // BOOST_VARIANT_DETAIL_SUBSTITUTE_FWD_HPP

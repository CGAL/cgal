//-----------------------------------------------------------------------------
// boost mpl/switch.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_SWITCH_HPP_INCLUDED
#define BOOST_MPL_SWITCH_HPP_INCLUDED

#include "boost/mpl/find_if.hpp"
#include "boost/mpl/select1st.hpp"
#include "boost/mpl/select2nd.hpp"
#include "boost/mpl/deref.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

namespace boost {
namespace mpl {

template< 
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Body)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct switch_
{
    typedef typename find_if<
          Body
        , apply1< lambda< select1st<_1> >, T >
        >::type iter_;
        
    typedef typename deref<iter_>::type pair_;
    typedef typename lambda< typename select2nd<pair_>::type >::type f_;
    typedef typename apply1<f_,T>::type type;

    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,switch_,(Body,T))
};

BOOST_MPL_AUX_VOID_SPEC(2, switch_)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_SWITCH_HPP_INCLUDED

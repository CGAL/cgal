
#ifndef BOOST_MPL_LAMBDA_FWD_HPP_INCLUDED
#define BOOST_MPL_LAMBDA_FWD_HPP_INCLUDED

// + file: boost/mpl/labmda_fwd.hpp
// + last modified: 02/aug/03

// Copyright (c) 2001-03
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/aux_/lambda_arity_param.hpp"
#include "boost/mpl/aux_/config/lambda.hpp"

namespace boost {
namespace mpl {

#if !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)

template< 
      typename T
    , typename Tag
    BOOST_MPL_AUX_LAMBDA_ARITY_PARAM(typename Arity)
    >
struct lambda;

#else

template< 
      typename T
    , typename Tag
    , bool Protect
    > 
struct lambda;

#endif

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_LAMBDA_FWD_HPP_INCLUDED

//-----------------------------------------------------------------------------
// boost mpl/always.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_ALWAYS_HPP_INCLUDED
#define BOOST_MPL_ALWAYS_HPP_INCLUDED

#include "boost/mpl/aux_/preprocessor/def_params_tail.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/arity_spec.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"

namespace boost {
namespace mpl {

template< typename Value >
struct always
{
    template<
          typename T
        BOOST_MPL_PP_NESTED_DEF_PARAMS_TAIL(1, typename T, void_)
        >
    struct apply
    {
        typedef Value type;
    };
};


BOOST_MPL_AUX_ARITY_SPEC(1,always)
BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(1,always)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_ALWAYS_HPP_INCLUDED

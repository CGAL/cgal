//-----------------------------------------------------------------------------
// boost mpl/aux_/empty_impl.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_EMPTY_IMPL_HPP_INCLUDED
#define BOOST_MPL_AUX_EMPTY_IMPL_HPP_INCLUDED

#include "boost/mpl/empty_fwd.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/aux_/traits_lambda_spec.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {

// default implementation; conrete sequences might override it by 
// specializing either the |empty_traits| or the primary |empty| template

template< typename Tag >
struct empty_traits
{
    template< typename Sequence > struct algorithm
        : is_same<
              typename begin<Sequence>::type
            , typename end<Sequence>::type
            >
    {
    };
};

BOOST_MPL_ALGORITM_TRAITS_LAMBDA_SPEC(1,empty_traits)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_EMPTY_IMPL_HPP_INCLUDED

//-----------------------------------------------------------------------------
// boost mpl/aux_/iter_apply.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_ITER_APPLY_HPP_INCLUDED
#define BOOST_MPL_ITER_APPLY_HPP_INCLUDED

#include "boost/mpl/apply.hpp"

namespace boost {
namespace mpl {
namespace aux {

template<
      typename F
    , typename Iterator
    >
struct iter_apply1
    : apply1<F,typename Iterator::type>
{
};

template<
      typename F
    , typename Iterator1
    , typename Iterator2
    >
struct iter_apply2
    : apply2<
          F
        , typename Iterator1::type
        , typename Iterator2::type
        >
{
};

} // namespace aux
} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_ITER_APPLY_HPP_INCLUDED

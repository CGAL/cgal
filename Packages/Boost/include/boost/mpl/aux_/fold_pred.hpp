//-----------------------------------------------------------------------------
// boost mpl/aux_/fold_pred.hpp header file
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

#ifndef BOOST_MPL_AUX_FOLD_PRED_HPP_INCLUDED
#define BOOST_MPL_AUX_FOLD_PRED_HPP_INCLUDED

#include "boost/mpl/same_as.hpp"
#include "boost/mpl/apply.hpp"

namespace boost {
namespace mpl {
namespace aux {

template< typename Last >
struct fold_pred
{
    template<
          typename State
        , typename Iterator
        >
    struct apply
        : not_same_as<Last>::template apply<Iterator>
    {
    };
};

} // namespace aux
} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_FOLD_PRED_HPP_INCLUDED

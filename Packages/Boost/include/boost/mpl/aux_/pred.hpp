//-----------------------------------------------------------------------------
// boost mpl/aux/pred.hpp header file
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

#ifndef BOOST_MPL_AUX_PRED_HPP_INCLUDED
#define BOOST_MPL_AUX_PRED_HPP_INCLUDED

namespace boost {
namespace mpl {
namespace aux {

// wrapper class to help users to deal with "legacy" metafunctions 
// (i.e. the ones that do not provide the '::type' interface)
//
// usage example: mpl::pred< boost::is_same<mpl::_1, int> >

template< typename Pred >
struct pred : Pred
{
    typedef pred type;
};

} // namespace aux
} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_PRED_HPP_INCLUDED

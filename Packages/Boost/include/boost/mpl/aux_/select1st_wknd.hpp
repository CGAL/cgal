//-----------------------------------------------------------------------------
// boost mpl/aux_/select1st_wknd.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_SELECT1ST_WKND_HPP_INCLUDED
#define BOOST_MPL_AUX_SELECT1ST_WKND_HPP_INCLUDED

#include "boost/mpl/aux_/is_msvc_eti_arg.hpp"
#include "boost/mpl/aux_/config/eti.hpp"

#if defined(BOOST_MPL_MSVC_ETI_BUG)

namespace boost { namespace mpl { namespace aux {

#   if defined(BOOST_MPL_MSVC_60_ETI_BUG)

template< typename Pair >
struct select1st_wknd
{
    typedef typename Pair::first type;
};

template<> struct select1st_wknd<int>
{
    typedef int type;
};

#   else // !defined(BOOST_MPL_MSVC_60_ETI_BUG)

template< bool > struct select1st_wknd_impl
{
    template< typename Pair > struct result_
    {
        typedef typename Pair::second type;
    };
};

template<> struct select1st_wknd_impl<false>
{
    template< typename Pair > struct result_
    {
        typedef int type;
    };
};

template< typename Pair > struct select1st_wknd
    : select1st_wknd_impl< !aux::is_msvc_eti_arg<Pair>::value >
        ::template result_<Pair>
{
};

#   endif // BOOST_MPL_MSVC_60_ETI_BUG workaround

}}} // namespace boost::mpl::aux

#   define BOOST_MPL_AUX_SELECT1ST_WKND(pair) ::boost::mpl::aux::select1st_wknd< pair >::type

#else // !defined(BOOST_MPL_MSVC_ETI_BUG)

#   define BOOST_MPL_AUX_SELECT1ST_WKND(pair) pair::first

#endif // BOOST_MPL_MSVC_ETI_BUG workaround

#endif // BOOST_MPL_AUX_SELECT1ST_WKND_HPP_INCLUDED

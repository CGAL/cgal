//-----------------------------------------------------------------------------
// boost mpl/aux_/msvc_msvc_eti_base.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
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

#ifndef BOOST_MPL_AUX_MSVC_ETI_BASE_HPP_INCLUDED
#define BOOST_MPL_AUX_MSVC_ETI_BASE_HPP_INCLUDED

#include "boost/mpl/aux_/config/eti.hpp"
#include "boost/mpl/aux_/is_msvc_eti_arg.hpp"

namespace boost { namespace mpl { namespace aux {

#if defined(BOOST_MPL_MSVC_ETI_BUG)

template< bool > struct msvc_eti_base_impl
{
    template< typename T > struct result_
    {
        typedef T type;
    };
};

template<> struct msvc_eti_base_impl<true>
{
    template< typename T > struct result_
    {
        typedef result_ type;
    };
};

template< typename T > struct msvc_eti_base
    : msvc_eti_base_impl< is_msvc_eti_arg<T>::value >
        ::template result_<T>
{
};

#else

template< typename T > struct msvc_eti_base
{
    typedef T type;
};

#endif // BOOST_MPL_MSVC_ETI_BUG

}}} // namespace boost::mpl::aux

#endif // BOOST_MPL_AUX_MSVC_ETI_BASE_HPP_INCLUDED

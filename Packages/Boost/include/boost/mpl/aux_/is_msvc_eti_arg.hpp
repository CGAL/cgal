//-----------------------------------------------------------------------------
// boost mpl/aux_/is_msvc_eti_arg.hpp header file
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

#ifndef BOOST_MPL_AUX_IS_MSVC_ETI_ARG_HPP_INCLUDED
#define BOOST_MPL_AUX_IS_MSVC_ETI_ARG_HPP_INCLUDED

#include "boost/mpl/aux_/yes_no.hpp"
#include "boost/mpl/aux_/config/eti.hpp"
#include "boost/mpl/aux_/config/static_constant.hpp"

namespace boost { namespace mpl { namespace aux {

#if defined(BOOST_MPL_MSVC_ETI_BUG)

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)

template< typename T >
struct is_msvc_eti_arg
{ 
    BOOST_STATIC_CONSTANT(bool, value = false);
};

#else

struct eti_int_convertible
{
    eti_int_convertible(int);
};

template< typename T >
struct is_msvc_eti_arg
{ 
    static no_tag test(...);
    static yes_tag test(eti_int_convertible);
    static T& get();

    BOOST_STATIC_CONSTANT(bool, value = 
          sizeof(test(get())) == sizeof(yes_tag)
        );
};

#endif // BOOST_MPL_MSVC_60_ETI_BUG

template<>
struct is_msvc_eti_arg<int>
{ 
    BOOST_STATIC_CONSTANT(bool, value = true);
};

#endif // BOOST_MPL_MSVC_ETI_BUG

}}} // namespace boost::mpl::aux

#endif // BOOST_MPL_AUX_IS_MSVC_ETI_ARG_HPP_INCLUDED

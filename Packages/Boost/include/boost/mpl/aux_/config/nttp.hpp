//-----------------------------------------------------------------------------
// boost mpl/aux_/config/nttp.hpp header file
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

#ifndef BOOST_MPL_AUX_CONFIG_NTTP_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_NTTP_HPP_INCLUDED

#include "boost/mpl/aux_/config/msvc.hpp"

// MSVC 6.5 ICE-s on the code as simple as this:
//
//    namespace std {
//    template< typename Char > struct string;
//    }
//
//    void foo(std::string<char>);
//
//    namespace boost { namespace mpl {
//    template< int > struct arg;
//    }}
//
// fortunately, a workaround is simple as well:
//
//    typedef int nttp_int;
//    template< nttp_int > struct arg;

#if defined(BOOST_MSVC) && BOOST_MSVC < 1300

#include "boost/preprocessor/cat.hpp"

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
namespace boost { namespace mpl {
typedef int     nttp_int;
typedef long    nttp_long;
}}
#endif

#   define BOOST_MPL_AUX_NTTP_DECL(T, x) BOOST_PP_CAT(nttp_,T) x /**/

#else
#   define BOOST_MPL_AUX_NTTP_DECL(T, x) T x /**/
#endif

#endif // BOOST_MPL_AUX_CONFIG_NTTP_HPP_INCLUDED

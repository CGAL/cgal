//-----------------------------------------------------------------------------
// boost mpl/aux_/msvc_dtw.hpp header file
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

// no include guards, the header is intended for multiple inclusion!

#include "boost/mpl/aux_/preprocessor/params.hpp"

// local macros, #undef-ined at the end of the header
#define AUX_DTW_PARAMS(param) \
    BOOST_MPL_PP_PARAMS(BOOST_MPL_AUX_MSVC_DTW_ARITY, param) \
/**/

#define AUX_DTW_ORIGINAL_NAME \
    BOOST_MPL_AUX_MSVC_DTW_ORIGINAL_NAME \
/**/

// warning: not a well-formed C++
// workaround for MSVC 6.5's "dependent template typedef bug"

template< typename F>
struct BOOST_MPL_AUX_MSVC_DTW_NAME
{
    template< bool > struct f_ : F {};
    template<> struct f_<true>
    {
        template< AUX_DTW_PARAMS(typename P) > struct AUX_DTW_ORIGINAL_NAME
        {
        };
    };

    template< AUX_DTW_PARAMS(typename T) > struct result_
        : f_< aux::msvc_never_true<F>::value >
            ::template AUX_DTW_ORIGINAL_NAME< AUX_DTW_PARAMS(T) >
    {
    };
};

#undef AUX_DTW_ORIGINAL_NAME
#undef AUX_DTW_PARAMS

#undef BOOST_MPL_AUX_MSVC_DTW_NAME
#undef BOOST_MPL_AUX_MSVC_DTW_ORIGINAL_NAME
#undef BOOST_MPL_AUX_MSVC_DTW_ARITY

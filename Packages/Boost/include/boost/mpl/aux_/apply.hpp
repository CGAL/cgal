//-----------------------------------------------------------------------------
// boost mpl/aux_/apply.hpp header file
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

#ifndef BOOST_MPL_AUX_APPLY_HPP_INCLUDED
#define BOOST_MPL_AUX_APPLY_HPP_INCLUDED

#include "boost/mpl/aux_/config/dtp.hpp"
#include "boost/config.hpp"

#define BOOST_MPL_AUX_APPLY(arity, args) \
    BOOST_PP_CAT(BOOST_MPL_AUX_APPLY,arity) args \
/**/

// agurt, 14/nov/02: temporary fix, need to research a couple of ICEs to
// get rid of this mess
#if defined(BOOST_MSVC) && BOOST_MSVC == 1300 && !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/apply.hpp"
#endif

#if defined(BOOST_MPL_USE_APPLY_INTERNALLY) \
    || defined(BOOST_BROKEN_DEFAULT_TEMPLATE_PARAMETERS_IN_NESTED_TEMPLATES) \
    || defined(BOOST_MSVC) && (BOOST_MSVC < 1300 || BOOST_MSVC == 1300 && defined(BOOST_MPL_PREPROCESSING_MODE))

#   if !defined(BOOST_MPL_PREPROCESSING_MODE)
#       include "boost/mpl/apply.hpp"
#   endif

// tokenization takes place before macro expansion (see 2.1 [lex.phases] 
// para 3-4), so, strictly speaking, spaces between '<', 'f', and '>' tokens 
// below (BOOST_MPL_AUX_APPLY0) are not required; they are needed in practice, 
// though, because there is at least one compiler (MSVC 6.5) that does not 
// conform to the standard here
#   define BOOST_MPL_AUX_APPLY0(f) apply0< f >
#   define BOOST_MPL_AUX_APPLY1(f,a1) apply1<f,a1>
#   define BOOST_MPL_AUX_APPLY2(f,a1,a2) apply2<f,a1,a2>
#   define BOOST_MPL_AUX_APPLY3(f,a1,a2,a3) apply3<f,a1,a2,a3>
#   define BOOST_MPL_AUX_APPLY4(f,a1,a2,a3,a4) apply4<f,a1,a2,a3,a4>
#   define BOOST_MPL_AUX_APPLY5(f,a1,a2,a3,a4,a5) apply5<f,a1,a2,a3,a4,a5>
#   define BOOST_MPL_AUX_APPLY6(f,a1,a2,a3,a4,a5,a6) apply6<f,a1,a2,a3,a4,a5,a6>
#   define BOOST_MPL_AUX_APPLY7(f,a1,a2,a3,a4,a5,a6,a7) apply7<f,a1,a2,a3,a4,a5,a6,a7>
#   define BOOST_MPL_AUX_APPLY8(f,a1,a2,a3,a4,a5,a6,a7,a8) apply8<f,a1,a2,a3,a4,a5,a6,a7,a8>
#   define BOOST_MPL_AUX_APPLY9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9) apply9<f,a1,a2,a3,a4,a5,a6,a7,a8,a9>

#else

#   define BOOST_MPL_AUX_APPLY0(f) f
#   define BOOST_MPL_AUX_APPLY1(f,a1) f::template apply<a1>
#   define BOOST_MPL_AUX_APPLY2(f,a1,a2) f::template apply<a1,a2>
#   define BOOST_MPL_AUX_APPLY3(f,a1,a2,a3) f::template apply<a1,a2,a3>
#   define BOOST_MPL_AUX_APPLY4(f,a1,a2,a3,a4) f::template apply<a1,a2,a3,a4>
#   define BOOST_MPL_AUX_APPLY5(f,a1,a2,a3,a4,a5) f::template apply<a1,a2,a3,a4,a5>
#   define BOOST_MPL_AUX_APPLY6(f,a1,a2,a3,a4,a5,a6) f::template apply<a1,a2,a3,a4,a5,a6>
#   define BOOST_MPL_AUX_APPLY7(f,a1,a2,a3,a4,a5,a6,a7) f::template apply<a1,a2,a3,a4,a5,a6,a7>
#   define BOOST_MPL_AUX_APPLY8(f,a1,a2,a3,a4,a5,a6,a7,a8) f::template apply<a1,a2,a3,a4,a5,a6,a7,a8>
#   define BOOST_MPL_AUX_APPLY9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9) f::template apply<a1,a2,a3,a4,a5,a6,a7,a8,a9>

#endif

#endif // BOOST_MPL_AUX_APPLY_HPP_INCLUDED

//-----------------------------------------------------------------------------
// boost mpl/aux_/arity_spec.hpp header file
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

#ifndef BOOST_MPL_AUX_ARITY_SPEC_HPP_INCLUDED
#define BOOST_MPL_AUX_ARITY_SPEC_HPP_INCLUDED

#include "boost/mpl/aux_/config/dtp.hpp"
#include "boost/mpl/aux_/preprocessor/params.hpp"
#include "boost/mpl/aux_/arity.hpp"
#include "boost/mpl/limits/arity.hpp"
#include "boost/config.hpp"

#if defined(BOOST_BROKEN_DEFAULT_TEMPLATE_PARAMETERS_IN_NESTED_TEMPLATES)
#   define BOOST_MPL_AUX_NONTYPE_ARITY_SPEC(i,type,name) \
namespace aux { \
template< BOOST_MPL_AUX_NTTP_DECL(int, N), BOOST_MPL_PP_PARAMS(i,type T) > \
struct arity< \
      name< BOOST_MPL_PP_PARAMS(i,T) > \
    , N \
    > \
{ \
    BOOST_STATIC_CONSTANT(int \
        , value = BOOST_MPL_METAFUNCTION_MAX_ARITY \
        ); \
}; \
} \
/**/
#else
#   define BOOST_MPL_AUX_NONTYPE_ARITY_SPEC(i,type,name) /**/
#endif

#   define BOOST_MPL_AUX_ARITY_SPEC(i,name) \
    BOOST_MPL_AUX_NONTYPE_ARITY_SPEC(i,typename,name) \
/**/

#endif // BOOST_MPL_AUX_ARITY_SPEC_HPP_INCLUDED

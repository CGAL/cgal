//-----------------------------------------------------------------------------
// boost mpl/aux_/void_spec.hpp header file
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

#ifndef BOOST_MPL_AUX_VOID_SPEC_HPP_INCLUDED
#define BOOST_MPL_AUX_VOID_SPEC_HPP_INCLUDED

#include "boost/mpl/lambda_fwd.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/int_fwd.hpp"
#include "boost/mpl/aux_/preprocessor/params.hpp"
#include "boost/mpl/aux_/preprocessor/enum.hpp"
#include "boost/mpl/aux_/preprocessor/def_params_tail.hpp"
#include "boost/mpl/aux_/arity.hpp"
#include "boost/mpl/aux_/template_arity_fwd.hpp"
#include "boost/mpl/aux_/lambda_arity_param.hpp"
#include "boost/mpl/aux_/algorithm_namespace.hpp"
#include "boost/mpl/aux_/config/dtp.hpp"
#include "boost/mpl/aux_/config/nttp.hpp"
#include "boost/mpl/aux_/config/ttp.hpp"
#include "boost/mpl/aux_/config/lambda.hpp"
#include "boost/mpl/aux_/config/overload_resolution.hpp"

#include "boost/config.hpp"

#define BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) \
    BOOST_MPL_PP_ENUM(i, void_) \
/**/

#if defined(BOOST_BROKEN_DEFAULT_TEMPLATE_PARAMETERS_IN_NESTED_TEMPLATES)
#   define BOOST_MPL_AUX_VOID_SPEC_ARITY(i, name) \
namespace aux { \
template< BOOST_MPL_AUX_NTTP_DECL(int, N) > \
struct arity< \
      name< BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) > \
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
#   define BOOST_MPL_AUX_VOID_SPEC_ARITY(i, name) /**/
#endif

#define BOOST_MPL_AUX_VOID_SPEC_MAIN(i, name) \
template<> \
struct name< BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) > \
{ \
    template< \
          BOOST_MPL_PP_PARAMS(i, typename T) \
        BOOST_MPL_PP_NESTED_DEF_PARAMS_TAIL(i, typename T, void_) \
        > \
    struct apply \
        : name< BOOST_MPL_PP_PARAMS(i, T) > \
    { \
    }; \
}; \
/**/

#if defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
#   define BOOST_MPL_AUX_VOID_SPEC_LAMBDA(i, name) \
template<> \
struct lambda< \
      name< BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) > \
    , void_ \
    , true \
    > \
{ \
    typedef name< BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) > type; \
}; \
/**/
#else
#   define BOOST_MPL_AUX_VOID_SPEC_LAMBDA(i, name) \
template<> \
struct lambda< \
      name< BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) > \
    , void_ \
    BOOST_MPL_AUX_LAMBDA_ARITY_PARAM(int_<-1>) \
    > \
{ \
    typedef name< BOOST_MPL_AUX_VOID_SPEC_PARAMS(i) > type; \
}; \
/**/
#endif

#if defined(BOOST_EXTENDED_TEMPLATE_PARAMETERS_MATCHING) || \
    defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT) && \
    defined(BOOST_MPL_BROKEN_OVERLOAD_RESOLUTION)
#   define BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(i, j, name) \
namespace aux { \
template< BOOST_MPL_PP_PARAMS(j, typename T) > \
struct template_arity< \
      name< BOOST_MPL_PP_PARAMS(j, T) > \
    > \
{ \
    BOOST_STATIC_CONSTANT(int, value = j); \
}; \
\
template<> \
struct template_arity< \
      name< BOOST_MPL_PP_ENUM(i, void_) > \
    > \
{ \
    BOOST_STATIC_CONSTANT(int, value = -1); \
}; \
} \
/**/
#else
#   define BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(i, j, name) /**/
#endif


#define BOOST_MPL_AUX_VOID_SPEC_PARAM(param) param = void_

#define BOOST_MPL_AUX_VOID_SPEC(i, name) \
BOOST_MPL_AUX_VOID_SPEC_MAIN(i, name) \
BOOST_MPL_AUX_VOID_SPEC_LAMBDA(i, name) \
BOOST_MPL_AUX_VOID_SPEC_ARITY(i, name) \
BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(i, i, name) \
/**/

#define BOOST_MPL_AUX_VOID_SPEC_EXT(i, j, name) \
BOOST_MPL_AUX_VOID_SPEC_MAIN(i, name) \
BOOST_MPL_AUX_VOID_SPEC_LAMBDA(i, name) \
BOOST_MPL_AUX_VOID_SPEC_ARITY(i, name) \
BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(i, j, name) \
/**/

#define BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(i, name) \
BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN \
BOOST_MPL_AUX_VOID_SPEC_MAIN(i, name) \
BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END \
BOOST_MPL_AUX_VOID_SPEC_LAMBDA(i, BOOST_MPL_AUX_AGLORITHM_NAMESPACE_PREFIX name) \
BOOST_MPL_AUX_VOID_SPEC_ARITY(i, BOOST_MPL_AUX_AGLORITHM_NAMESPACE_PREFIX name) \
BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(i, i, BOOST_MPL_AUX_AGLORITHM_NAMESPACE_PREFIX name) \
/**/

#endif // BOOST_MPL_AUX_VOID_SPEC_HPP_INCLUDED

//-----------------------------------------------------------------------------
// boost mpl/apply.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-03
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_APPLY_HPP_INCLUDED
#define BOOST_MPL_APPLY_HPP_INCLUDED

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/arg_fwd.hpp"
#   include "boost/mpl/void.hpp"
#   include "boost/mpl/aux_/arity.hpp"
#   include "boost/mpl/aux_/msvc_never_true.hpp"
#   include "boost/type_traits/same_traits.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
 && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER apply.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/limits/arity.hpp"
#   include "boost/mpl/aux_/lambda_support.hpp"
#   include "boost/mpl/aux_/preprocessor/params.hpp"
#   include "boost/mpl/aux_/preprocessor/default_params.hpp"
#   include "boost/mpl/aux_/preprocessor/partial_spec_params.hpp"
#   include "boost/mpl/aux_/preprocessor/enum.hpp"
#   include "boost/mpl/aux_/preprocessor/add.hpp"
#   include "boost/mpl/aux_/config/dtp.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"
#   include "boost/mpl/aux_/config/eti.hpp"
#   include "boost/mpl/aux_/config/lambda.hpp"

#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/iterate.hpp"

#   include "boost/config.hpp"

// agurt, 15/jan/02: top-level 'apply' template gives an ICE on MSVC
// (for known reasons)
#if defined(BOOST_MSVC) && (BOOST_MSVC <= 1200)
#   define BOOST_MPL_NO_APPLY_TEMPLATE
#endif

namespace boost {
namespace mpl {

// local macros, #undef-ined at the end of the header
#   define AUX_APPLY_PARAMS(param) \
    BOOST_MPL_PP_PARAMS( \
          BOOST_MPL_METAFUNCTION_MAX_ARITY \
        , param \
        ) \
    /**/

#   define AUX_APPLY_DEFAULT_PARAMS(param, value) \
    BOOST_MPL_PP_DEFAULT_PARAMS( \
          BOOST_MPL_METAFUNCTION_MAX_ARITY \
        , param \
        , value \
        ) \
    /**/

#   define AUX_APPLY_N_PARAMS(n, param) \
    BOOST_MPL_PP_PARAMS(n, param) \
    /**/

#   define AUX_APPLY_N_COMMA_PARAMS(n, param) \
    BOOST_PP_COMMA_IF(n) \
    BOOST_MPL_PP_PARAMS(n, param) \
    /**/

#   define AUX_APPLY_N_PARTIAL_SPEC_PARAMS(n, param, def) \
    BOOST_PP_COMMA_IF(n) \
    BOOST_MPL_PP_PARTIAL_SPEC_PARAMS(n, param, def) \
    /**/
    
#   define AUX_APPLY_N_SPEC_PARAMS(n, param) \
    BOOST_MPL_PP_ENUM(BOOST_PP_INC(n), param) \
    /**/

#   if !defined(BOOST_MPL_NO_APPLY_TEMPLATE)

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
// forward declaration
template<
      typename F, AUX_APPLY_DEFAULT_PARAMS(typename T, void_)
    >
struct apply;
#else
namespace aux {
template< BOOST_MPL_AUX_NTTP_DECL(int, arity_) > struct apply_impl_chooser;
}
#endif

#   endif // BOOST_MPL_NO_APPLY_TEMPLATE

#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_METAFUNCTION_MAX_ARITY, "boost/mpl/apply.hpp"))
#include BOOST_PP_ITERATE()

#   if !defined(BOOST_MPL_NO_APPLY_TEMPLATE)
// real C++ version is already taken care of
#   if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace aux {
// apply_count_args
#define BOOST_MPL_AUX_COUNT_ARGS_PREFIX apply
#define BOOST_MPL_AUX_COUNT_ARGS_DEFAULT void_
#define BOOST_MPL_AUX_COUNT_ARGS_ARITY BOOST_MPL_METAFUNCTION_MAX_ARITY
#include "boost/mpl/aux_/count_args.hpp"
} // namespace aux

template<
      typename F, AUX_APPLY_DEFAULT_PARAMS(typename T, void_)
    >
struct apply
    : aux::apply_impl_chooser< 
          aux::apply_count_args< AUX_APPLY_PARAMS(T) >::value
        >::template result_< F, AUX_APPLY_PARAMS(T) >::type
{
};

#   endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
#   endif // BOOST_MPL_NO_APPLY_TEMPLATE

#   undef AUX_APPLY_N_SPEC_PARAMS
#   undef AUX_APPLY_N_PARTIAL_SPEC_PARAMS
#   undef AUX_APPLY_N_COMMA_PARAMS
#   undef AUX_APPLY_N_PARAMS
#   undef AUX_APPLY_DEFAULT_PARAMS
#   undef AUX_APPLY_PARAMS

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_APPLY_HPP_INCLUDED

///// iteration, depth == 1

#elif BOOST_PP_ITERATION_DEPTH() == 1

#   define i BOOST_PP_FRAME_ITERATION(1)
#   if i == 0

template< typename F >
struct apply0 : F
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1, apply0, (F))
};

#if defined(BOOST_MPL_MSVC_ETI_BUG)
//: workaround for the ETI bug
template<>
struct apply0<int>
{
    typedef int type;
};
#endif

#   else // i > 0

#   if defined(BOOST_MSVC) && (BOOST_MSVC < 1300)
// MSVC version

namespace aux {
// msvc_apply##i
#define BOOST_MPL_AUX_MSVC_DTW_NAME BOOST_PP_CAT(msvc_apply,i)
#define BOOST_MPL_AUX_MSVC_DTW_ORIGINAL_NAME apply
#define BOOST_MPL_AUX_MSVC_DTW_ARITY i
#include "boost/mpl/aux_/msvc_dtw.hpp"
} // namespace aux

template<
      typename F, AUX_APPLY_N_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(apply,i)
{
    // Metafunction forwarding confuses vc6
    typedef typename BOOST_PP_CAT(aux::msvc_apply,i)<F>::template result_<
        AUX_APPLY_N_PARAMS(i, T)
    >::type type;
    
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          BOOST_PP_INC(i)
        , BOOST_PP_CAT(apply,i)
        , (F, AUX_APPLY_N_PARAMS(i,T))
        )
};

#   elif defined(BOOST_BROKEN_DEFAULT_TEMPLATE_PARAMETERS_IN_NESTED_TEMPLATES)
// MWCW/Borland version

namespace aux {
template<
      int N, typename F, AUX_APPLY_N_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(apply_impl,i);
}

#define BOOST_PP_ITERATION_PARAMS_2 \
    (3,(0, BOOST_MPL_METAFUNCTION_MAX_ARITY - i, "boost/mpl/apply.hpp"))
#include BOOST_PP_ITERATE()

template<
      typename F, AUX_APPLY_N_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(apply,i)
    : BOOST_PP_CAT(aux::apply_impl,i)<
          ::boost::mpl::aux::arity<F,i>::value
        , F
        , AUX_APPLY_N_PARAMS(i, T)
        >::type
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          BOOST_PP_INC(i)
        , BOOST_PP_CAT(apply,i)
        , (F, AUX_APPLY_N_PARAMS(i,T))
        )
};

#   else
// ISO98 C++, with minor concession to vc7

template<
      typename F, AUX_APPLY_N_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(apply,i)
{
    // Metafunction forwarding confuses vc7
    typedef typename F::template apply<
       AUX_APPLY_N_PARAMS(i, T)
    >::type type;
        
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          BOOST_PP_INC(i)
        , BOOST_PP_CAT(apply,i)
        , (F, AUX_APPLY_N_PARAMS(i,T))
        )
};

#   endif // workarounds

#if defined(BOOST_MPL_MSVC_ETI_BUG)
//: workaround for ETI bug
template<>
struct BOOST_PP_CAT(apply,i)<AUX_APPLY_N_SPEC_PARAMS(i, int)>
{
    typedef int type;
};
#endif

#   endif // i > 0

#   if !defined(BOOST_MPL_NO_APPLY_TEMPLATE)
#   if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#if i == BOOST_MPL_METAFUNCTION_MAX_ARITY

//: primary template (not a specialization!)
template<
      typename F AUX_APPLY_N_COMMA_PARAMS(i, typename T)
    >
struct apply
    : BOOST_PP_CAT(apply,i)< F AUX_APPLY_N_COMMA_PARAMS(i, T) >
{
};

#else

template<
      typename F AUX_APPLY_N_COMMA_PARAMS(i, typename T)
    >
struct apply< F AUX_APPLY_N_PARTIAL_SPEC_PARAMS(i, T, void_) >
    : BOOST_PP_CAT(apply,i)< F AUX_APPLY_N_COMMA_PARAMS(i, T) >
{
};

#endif // i == BOOST_MPL_METAFUNCTION_MAX_ARITY

#   else

namespace aux {

template<>
struct apply_impl_chooser<i>
{
    template<
          typename F, AUX_APPLY_PARAMS(typename T)
        >
    struct result_
    {
        typedef BOOST_PP_CAT(apply,i)<
              F AUX_APPLY_N_COMMA_PARAMS(i, T)
            > type;
    };
};

} // namespace aux

#   endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
#   endif // BOOST_MPL_NO_APPLY_TEMPLATE

#   undef i

///// iteration, depth == 2

#elif BOOST_PP_ITERATION_DEPTH() == 2

#   define j BOOST_PP_FRAME_ITERATION(2)

namespace aux {

template<
      typename F, AUX_APPLY_N_PARAMS(i, typename T)
    >
struct BOOST_PP_CAT(apply_impl,i)<
          BOOST_MPL_PP_ADD(i, j)
        , F
        , AUX_APPLY_N_PARAMS(i, T)
        >
{
    typedef typename F::template apply<
          AUX_APPLY_N_PARAMS(i, T)
        BOOST_PP_COMMA_IF(j) BOOST_MPL_PP_ENUM(j, void_)
        > type;
};

} // namespace aux

#   undef j

#endif // BOOST_PP_IS_ITERATING

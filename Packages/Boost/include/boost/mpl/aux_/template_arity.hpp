//-----------------------------------------------------------------------------
// boost mpl/aux_/template_arity.hpp header file
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

#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_AUX_TEMPLATE_ARITY_HPP_INCLUDED
#define BOOST_MPL_AUX_TEMPLATE_ARITY_HPP_INCLUDED

#include "boost/mpl/aux_/config/ttp.hpp"
#include "boost/mpl/aux_/config/lambda.hpp"

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/aux_/template_arity_fwd.hpp"
#   if !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
#   if defined(BOOST_EXTENDED_TEMPLATE_PARAMETERS_MATCHING)
#       include "boost/mpl/aux_/type_wrapper.hpp"
#   endif
#   else
#       include "boost/mpl/aux_/has_rebind.hpp"
#   endif
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
 && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER template_arity.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   if !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
#   if defined(BOOST_EXTENDED_TEMPLATE_PARAMETERS_MATCHING)

#   include "boost/mpl/limits/arity.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"
#   include "boost/mpl/aux_/preprocessor/range.hpp"
#   include "boost/mpl/aux_/preprocessor/repeat.hpp"
#   include "boost/mpl/aux_/preprocessor/params.hpp"

#   include "boost/preprocessor/seq/fold_left.hpp"
#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

namespace boost { namespace mpl { namespace aux {

template< BOOST_MPL_AUX_NTTP_DECL(int, N) > struct arity_tag
{
    typedef char (&type)[N + 1];
};

#define AUX_MAX_ARITY_OP(unused, state, i) \
    ( BOOST_PP_CAT(C,i) > 0 ? BOOST_PP_CAT(C,i) : state ) \
/**/

template<
      BOOST_MPL_PP_PARAMS(
          BOOST_MPL_METAFUNCTION_MAX_ARITY
        , BOOST_MPL_AUX_NTTP_DECL(int, C)
        )
    >
struct max_arity
{
    BOOST_STATIC_CONSTANT(int, value = 
          BOOST_PP_SEQ_FOLD_LEFT(
              AUX_MAX_ARITY_OP
            , -1
            , BOOST_MPL_PP_RANGE(1, BOOST_MPL_METAFUNCTION_MAX_ARITY)
            )
        );
};

#undef AUX_MAX_ARITY_OP

arity_tag<0> arity_helper(...);

#define BOOST_PP_ITERATION_LIMITS (1, BOOST_MPL_METAFUNCTION_MAX_ARITY)
#define BOOST_PP_FILENAME_1 "boost/mpl/aux_/template_arity.hpp"
#include BOOST_PP_ITERATE()

template< typename F, BOOST_MPL_AUX_NTTP_DECL(int, N) >
struct template_arity_impl
{
    BOOST_STATIC_CONSTANT(int, value = 
          sizeof(arity_helper(type_wrapper<F>(),arity_tag<N>())) - 1
        );
};

#define AUX_TEMPLATE_ARITY_IMPL_INVOCATION(unused, i, F) \
    BOOST_PP_COMMA_IF(i) template_arity_impl<F,BOOST_PP_INC(i)>::value \
/**/

template< typename F >
struct template_arity
{
    BOOST_STATIC_CONSTANT(int, value = (
          max_arity< BOOST_MPL_PP_REPEAT(
              BOOST_MPL_METAFUNCTION_MAX_ARITY
            , AUX_TEMPLATE_ARITY_IMPL_INVOCATION
            , F
            ) >::value
        ));
};

#undef AUX_TEMPLATE_ARITY_IMPL_INVOCATION

}}} // namespace boost::mpl::aux

#   endif // BOOST_EXTENDED_TEMPLATE_PARAMETERS_MATCHING
#   else // BOOST_MPL_NO_FULL_LAMBDA_SUPPORT

#   include "boost/mpl/aux_/config/eti.hpp"
#   include "boost/mpl/aux_/config/static_constant.hpp"
#   include "boost/mpl/aux_/config/workaround.hpp"

namespace boost { namespace mpl { namespace aux {

template< bool >
struct template_arity_impl
{
    template< typename F > struct result_
    {
        BOOST_STATIC_CONSTANT(int, value = -1);
    };
};

template<>
struct template_arity_impl<true>
{
    template< typename F > struct result_
    {
#if defined(__BORLANDC__) && (__BORLANDC__ >= 0x561 && !defined(BOOST_STRICT_CONFIG))
        enum { value = F::arity };
#else
        BOOST_STATIC_CONSTANT(int, value = F::arity);
#endif
    };
};

template< typename F >
struct template_arity
    : template_arity_impl< ::boost::mpl::aux::has_rebind<F>::value >
        ::template result_<F>
{
};

#if defined(BOOST_MPL_MSVC_ETI_BUG)
template<>
struct template_arity<int>
{
    BOOST_STATIC_CONSTANT(int, value = -1);
};
#endif

}}} // namespace boost::mpl::aux

#   endif // BOOST_MPL_NO_FULL_LAMBDA_SUPPORT

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_TEMPLATE_ARITY_HPP_INCLUDED

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

template<
      template< BOOST_MPL_PP_PARAMS(i, typename P) > class F
    , BOOST_MPL_PP_PARAMS(i, typename T)
    >
typename arity_tag<i>::type
arity_helper(type_wrapper< F<BOOST_MPL_PP_PARAMS(i, T)> >, arity_tag<i>);

#undef i
#endif // BOOST_PP_IS_ITERATING

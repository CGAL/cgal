//-----------------------------------------------------------------------------
// boost mpl/inherit.hpp header file
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

#ifndef BOOST_MPL_INHERIT_HPP_INCLUDED
#define BOOST_MPL_INHERIT_HPP_INCLUDED

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/empty_base.hpp"
#   include "boost/mpl/aux_/void_spec.hpp"
#   include "boost/mpl/aux_/lambda_support.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if 0 //!defined(BOOST_MPL_NO_PREPROCESSED_HEADERS)
 //&& !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER inherit.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/aux_/preprocessor/params.hpp"
#   include "boost/mpl/aux_/preprocessor/default_params.hpp"
#   include "boost/mpl/limits/arity.hpp"

#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/dec.hpp"
#   include "boost/preprocessor/cat.hpp"
#   include "boost/config.hpp"

namespace boost {
namespace mpl {

// 'inherit<T1,T2,..,Tn>' metafunction; returns an unspecified class type
// produced by public derivation from all metafunction's parameters 
// (T1,T2,..,Tn), except the parameters of 'empty_base' class type; 
// regardless the position and number of 'empty_base' parameters in the 
// metafunction's argument list, derivation from them is always a no-op;
// for instance:
//      inherit<her>::type == her
//      inherit<her,my>::type == struct unspecified : her, my {};
//      inherit<empty_base,her>::type == her
//      inherit<empty_base,her,empty_base,empty_base>::type == her
//      inherit<her,empty_base,my>::type == struct unspecified : her, my {};
//      inherit<empty_base,empty_base>::type == empty_base

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template< 
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    > 
struct inherit2
    : T1, T2
{
    typedef inherit2 type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2, inherit2, (T1,T2))
};

template< typename T1 > 
struct inherit2<T1,empty_base>
{
    typedef T1 type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(2, inherit2, (T1,empty_base))
};

template< typename T2 > 
struct inherit2<empty_base,T2>
{
    typedef T2 type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(2, inherit2, (empty_base,T2))
};

// needed to disambiguate the previous two in case when both 
// T1 and T2 == empty_base
template<> 
struct inherit2<empty_base,empty_base>
{
    typedef empty_base type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(2, inherit2, (empty_base,empty_base))
};

#else

namespace aux {

template< bool C1, bool C2 >
struct inherit2_impl
{
    template< typename Derived, typename T1, typename T2 > struct result_ 
        : T1, T2
    {
        typedef Derived type_;
    };
};

template<>
struct inherit2_impl<false,true>
{
    template< typename Derived, typename T1, typename T2 > struct result_
        : T1
    {
        typedef T1 type_;
    };
};

template<>
struct inherit2_impl<true,false>
{
    template< typename Derived, typename T1, typename T2 > struct result_
        : T2 
    {
        typedef T2 type_;
    };
};

template<>
struct inherit2_impl<true,true>
{
    template< typename Derived, typename T1, typename T2 > struct result_
    {
        typedef T1 type_;
    };
};

} // namespace aux

template< 
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    > 
struct inherit2
    : aux::inherit2_impl<
          is_empty_base<T1>::value
        , is_empty_base<T2>::value
        >::template result_< inherit2<T1,T2>,T1,T2 >
{
    typedef typename inherit2::type_ type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2, inherit2, (T1,T2))
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

BOOST_MPL_AUX_VOID_SPEC(2, inherit2)

#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(3, BOOST_MPL_METAFUNCTION_MAX_ARITY, "boost/mpl/inherit.hpp"))
#include BOOST_PP_ITERATE()

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_INHERIT_HPP_INCLUDED

///// iteration

#else
#define n BOOST_PP_FRAME_ITERATION(1)

template<
      BOOST_MPL_PP_DEFAULT_PARAMS(n, typename T, void_)
    >
struct BOOST_PP_CAT(inherit,n)
    : inherit2<
          typename BOOST_PP_CAT(inherit,BOOST_PP_DEC(n))<
              BOOST_MPL_PP_PARAMS(BOOST_PP_DEC(n), T)
            >::type
        , BOOST_PP_CAT(T,n)
        >
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          n
        , BOOST_PP_CAT(inherit,n)
        , (BOOST_MPL_PP_PARAMS(n, T))
        )
};

BOOST_MPL_AUX_VOID_SPEC(n, BOOST_PP_CAT(inherit,n))

#if n == BOOST_MPL_METAFUNCTION_MAX_ARITY
//: primary template
template<
      BOOST_MPL_PP_DEFAULT_PARAMS(n, typename T, empty_base)
    >
struct inherit
    : BOOST_PP_CAT(inherit,n)<BOOST_MPL_PP_PARAMS(n, T)>
{
};

// 'void_' specialization
template<>
struct inherit< BOOST_MPL_AUX_VOID_SPEC_PARAMS(n) >
{
    template<
          BOOST_MPL_PP_NESTED_DEF_PARAMS_TAIL(0, typename T, empty_base)
        >
    struct apply
        : inherit< BOOST_MPL_PP_PARAMS(n, T) >
    {
    };
};

BOOST_MPL_AUX_VOID_SPEC_LAMBDA(n, inherit)
BOOST_MPL_AUX_VOID_SPEC_ARITY(n, inherit)
BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(n, n, inherit)
#endif

#undef n
#endif // BOOST_PP_IS_ITERATING

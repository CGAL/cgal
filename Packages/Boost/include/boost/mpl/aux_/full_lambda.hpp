
#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_AUX_FULL_LAMBDA_HPP_INCLUDED
#define BOOST_MPL_AUX_FULL_LAMBDA_HPP_INCLUDED

// + file: boost/mpl/aux_/full_lambda.hpp
// + last modified: 03/aug/03

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
//
// See http://www.boost.org/libs/mpl for documentation.

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/lambda_fwd.hpp"
#   include "boost/mpl/bind.hpp"
#   include "boost/mpl/protect.hpp"
#   include "boost/mpl/quote.hpp"
#   include "boost/mpl/bool.hpp"
#   include "boost/mpl/int_fwd.hpp"
#   include "boost/mpl/aux_/template_arity.hpp"
#   include "boost/mpl/aux_/config/ttp.hpp"
#endif

#include "boost/mpl/aux_/lambda_expr.hpp"
#include "boost/mpl/aux_/lambda_arity_param.hpp"
#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
 && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER full_lambda.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/limits/arity.hpp"
#   include "boost/mpl/aux_/preprocessor/default_params.hpp"
#   include "boost/mpl/aux_/preprocessor/params.hpp"
#   include "boost/mpl/aux_/preprocessor/enum.hpp"
#   include "boost/mpl/aux_/preprocessor/repeat.hpp"

#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

namespace boost {
namespace mpl {

// local macros, #undef-ined at the end of the header
#   define AUX_LAMBDA_PARAMS(i, param) \
    BOOST_MPL_PP_PARAMS(i, param) \
    /**/

#   define AUX_LAMBDA_BIND_PARAMS(param) \
    BOOST_MPL_PP_PARAMS( \
          BOOST_MPL_METAFUNCTION_MAX_ARITY \
        , param \
        ) \
    /**/

#   define AUX_LAMBDA_BIND_N_PARAMS(i, param) \
    BOOST_PP_COMMA_IF(i) \
    BOOST_MPL_PP_PARAMS(i, param) \
    /**/

#   define AUX_ARITY_PARAM(param) \
    BOOST_MPL_AUX_LAMBDA_ARITY_PARAM(param) \
    /**/

template<
      typename T
    , typename Tag
    , typename Protect = false_
    AUX_ARITY_PARAM(typename Arity = int_< aux::template_arity<T>::value >)
    >
struct lambda_impl
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(false_)
    typedef T type;
};

template<
      typename T
    , typename Tag = void_
    AUX_ARITY_PARAM(typename Arity = int_< aux::template_arity<T>::value >)
    >
struct lambda
    : lambda_impl<T,Tag,false_ AUX_ARITY_PARAM(Arity)>
{
};

#if !defined(BOOST_MPL_NO_LAMBDA_HEURISTIC)

#define n BOOST_MPL_METAFUNCTION_MAX_ARITY
namespace aux {

template<
      BOOST_MPL_PP_DEFAULT_PARAMS(n,bool C,false)
    >
struct lambda_or
    : true_
{
};

template<>
struct lambda_or< BOOST_MPL_PP_ENUM(n,false) >
    : false_
{
};

} // namespace aux
#undef n

template< int N, typename Tag, typename Protect >
struct lambda_impl< arg<N>,Tag,Protect AUX_ARITY_PARAM(int_<-1>) >
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(true_)
    typedef mpl::arg<N> type; // qualified for the sake of MIPSpro 7.41
};

#endif // BOOST_MPL_NO_LAMBDA_HEURISTIC

#define BOOST_PP_ITERATION_LIMITS (0, BOOST_MPL_METAFUNCTION_MAX_ARITY)
#define BOOST_PP_FILENAME_1 "boost/mpl/aux_/full_lambda.hpp"
#include BOOST_PP_ITERATE()

//: special case for 'protect'
template< typename T, typename Tag, typename Protect >
struct lambda_impl< protect<T>,Tag,Protect AUX_ARITY_PARAM(int_<1>) >
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(false_)
    typedef protect<T> type;
};

//: specializations for main 'bind', 'bind1st' and 'bind2nd' forms
template<
      typename F, AUX_LAMBDA_BIND_PARAMS(typename T)
    , typename Tag
    , typename Protect
    >
struct lambda_impl<
      bind<F,AUX_LAMBDA_BIND_PARAMS(T)>
    , Tag
    , Protect 
    AUX_ARITY_PARAM(int_<BOOST_PP_INC(BOOST_MPL_METAFUNCTION_MAX_ARITY)>)
    >
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(false_)
    typedef bind<F, AUX_LAMBDA_BIND_PARAMS(T)> type;
};

template<
      typename F, typename T
    , typename Tag
    , typename Protect
    >
struct lambda_impl< bind1st<F,T>,Tag,Protect AUX_ARITY_PARAM(int_<2>) >
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(false_)
    typedef bind1st<F,T> type;
};

template<
      typename F, typename T
    , typename Tag
    , typename Protect
    >
struct lambda_impl< bind2nd<F,T>,Tag,Protect AUX_ARITY_PARAM(int_<2>) >
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(false_)
    typedef bind2nd<F,T> type;
};

#   undef AUX_ARITY_PARAM
#   undef AUX_LAMBDA_BIND_N_PARAMS
#   undef AUX_LAMBDA_BIND_PARAMS
#   undef AUX_LAMBDA_PARAMS

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_FULL_LAMBDA_HPP_INCLUDED

///// iteration, depth == 1

#elif BOOST_PP_ITERATION_DEPTH() == 1
#define i BOOST_PP_FRAME_ITERATION(1)

#if i > 0

template<
      template< AUX_LAMBDA_PARAMS(i, typename P) > class F
    , AUX_LAMBDA_PARAMS(i, typename T)
    , typename Tag
    >
struct lambda< F<AUX_LAMBDA_PARAMS(i, T)>, Tag AUX_ARITY_PARAM(int_<i>) >
    : lambda_impl< F<AUX_LAMBDA_PARAMS(i, T)>,Tag,true_ AUX_ARITY_PARAM(int_<i>) >
{
};

#if defined(BOOST_MPL_NO_LAMBDA_HEURISTIC)

template<
      template< AUX_LAMBDA_PARAMS(i, typename P) > class F
    , AUX_LAMBDA_PARAMS(i, typename T)
    , typename Tag
    , typename Protect
    >
struct lambda_impl< 
      F<AUX_LAMBDA_PARAMS(i, T)>, Tag, Protect AUX_ARITY_PARAM(int_<i>)
    >
{
#   define AUX_LAMBDA_INVOCATION(unused, i, T) \
    BOOST_PP_COMMA_IF(i) \
    typename lambda_impl< BOOST_PP_CAT(T, BOOST_PP_INC(i)),Tag >::type \
    /**/

    typedef BOOST_PP_CAT(bind,i)<
          BOOST_PP_CAT(quote,i)<F,Tag>
        , BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_INVOCATION, T)
        > type;

#   undef AUX_LAMBDA_INVOCATION
};

#else // BOOST_MPL_NO_LAMBDA_HEURISTIC

#   define AUX_LAMBDA_RESULT(unused, i, T) \
    BOOST_PP_COMMA_IF(i) \
    typename BOOST_PP_CAT(T, BOOST_PP_INC(i))::type \
    /**/

namespace aux {

template<
      typename IsLE, typename Tag, typename Protect
    , template< AUX_LAMBDA_PARAMS(i, typename P) > class F
    , AUX_LAMBDA_PARAMS(i, typename L)
    >
struct BOOST_PP_CAT(le_result,i)
{
    typedef F<
          BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_RESULT, L)
        > type;
};

template<
      typename Tag
    , template< AUX_LAMBDA_PARAMS(i, typename P) > class F
    , AUX_LAMBDA_PARAMS(i, typename L)
    >
struct BOOST_PP_CAT(le_result,i)< true_,Tag,false_,F,AUX_LAMBDA_PARAMS(i, L) >
{
    typedef BOOST_PP_CAT(bind,i)<
          BOOST_PP_CAT(quote,i)<F,Tag>
        , BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_RESULT, L)
        > type;
};

template<
      typename Tag
    , template< AUX_LAMBDA_PARAMS(i, typename P) > class F
    , AUX_LAMBDA_PARAMS(i, typename L)
    >
struct BOOST_PP_CAT(le_result,i)< true_,Tag,true_,F,AUX_LAMBDA_PARAMS(i, L) >
{
    typedef protect< BOOST_PP_CAT(bind,i)<
          BOOST_PP_CAT(quote,i)<F,Tag>
        , BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_RESULT, L)
        > > type;
};

} // namespace aux

#   define AUX_LAMBDA_INVOCATION(unused, i, T) \
    typedef lambda_impl< BOOST_PP_CAT(T, BOOST_PP_INC(i)), Tag > \
        BOOST_PP_CAT(l,BOOST_PP_INC(i)); \
    /**/

#   define AUX_IS_LAMBDA_EXPR(unused, i, unused2) \
    BOOST_PP_COMMA_IF(i) \
    BOOST_PP_CAT(l,BOOST_PP_INC(i))::is_le::value \
    /**/

template<
      template< AUX_LAMBDA_PARAMS(i, typename P) > class F
    , AUX_LAMBDA_PARAMS(i, typename T)
    , typename Tag
    , typename Protect
    >
struct lambda_impl< 
      F<AUX_LAMBDA_PARAMS(i, T)>, Tag, Protect AUX_ARITY_PARAM(int_<i>)
    >
{
    BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_INVOCATION, T)
    typedef aux::lambda_or<
          BOOST_MPL_PP_REPEAT(i, AUX_IS_LAMBDA_EXPR, unused)
        > is_le;

    typedef typename aux::BOOST_PP_CAT(le_result,i)<
          typename is_le::type
        , Tag
        , Protect
        , F
        , AUX_LAMBDA_PARAMS(i, l)
        >::type type;
};

#   undef AUX_IS_LAMBDA_EXPR
#   undef AUX_LAMBDA_INVOCATION
#   undef AUX_LAMBDA_RESULT

#endif // BOOST_MPL_NO_LAMBDA_HEURISTIC
#endif // i > 0

template<
      typename F AUX_LAMBDA_BIND_N_PARAMS(i, typename T)
    , typename Tag
    , typename Protect
    >
struct lambda_impl<
      BOOST_PP_CAT(bind,i)<F AUX_LAMBDA_BIND_N_PARAMS(i, T)>
    , Tag
    , Protect AUX_ARITY_PARAM(int_<BOOST_PP_INC(i)>)
    >
{
    BOOST_MPL_AUX_IS_LAMBDA_EXPR(false_)
    typedef BOOST_PP_CAT(bind,i)<
          F
        AUX_LAMBDA_BIND_N_PARAMS(i, T)
        > type;
};

#undef i
#endif // BOOST_PP_IS_ITERATING


#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_AUX_LAMBDA_NO_CTPS_HPP_INCLUDED
#define BOOST_MPL_AUX_LAMBDA_NO_CTPS_HPP_INCLUDED

// + file: boost/mpl/aux_/lambda_no_ctps.hpp
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
#   include "boost/mpl/bool.hpp"
#   include "boost/mpl/void.hpp"
#   include "boost/mpl/aux_/template_arity.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if    !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER lambda_no_ctps.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/aux_/config/nttp.hpp"
#   include "boost/mpl/limits/arity.hpp"
#   include "boost/mpl/aux_/preprocessor/params.hpp"
#   include "boost/mpl/aux_/preprocessor/repeat.hpp"

#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

#   include "boost/config.hpp"

namespace boost {
namespace mpl {

// local macros, #undef-ined at the end of the header
#   define AUX_LAMBDA_PARAMS(i, param) \
    BOOST_MPL_PP_PARAMS(i, param) \
    /**/

namespace aux {

template< BOOST_MPL_AUX_NTTP_DECL(int, arity_), bool Protect > struct lambda_impl
{
    template< typename T, typename Tag > struct result_
    {
        typedef T type;
    };
};

#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(1, BOOST_MPL_METAFUNCTION_MAX_ARITY, "boost/mpl/aux_/lambda_no_ctps.hpp"))
#include BOOST_PP_ITERATE()

} // namespace aux

template<
      typename T
    , typename Tag = void_
    , bool Protect = true 
    >
struct lambda
{
    // Metafunction forwarding confuses vc6
    typedef typename aux::lambda_impl<
          ::boost::mpl::aux::template_arity<T>::value
#if !defined(BOOST_MSVC) || BOOST_MSVC > 1200
        , Protect
#else
        , bool_<Protect>::value
#endif
    >::template result_<T,Tag>::type type;
};

#   undef AUX_LAMBDA_PARAMS

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_LAMBDA_NO_CTPS_HPP_INCLUDED

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

#   define AUX_LAMBDA_INVOCATION(unused, i, T) \
    typedef typename lambda< \
          typename F::BOOST_PP_CAT(arg,BOOST_PP_INC(i)) \
        , Tag \
        , false \
        >::type BOOST_PP_CAT(arg,BOOST_PP_INC(i)); \
    /**/

template<> struct lambda_impl<i,false>
{
    template< typename F, typename Tag > struct result_
    {
        typedef typename F::rebind f_;
        BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_INVOCATION, T)

        typedef BOOST_PP_CAT(bind,i)<
              f_
            , AUX_LAMBDA_PARAMS(i, arg)
            > type;
    };
};

template<> struct lambda_impl<i,true>
{
    template< typename F, typename Tag > struct result_
    {
        typedef typename F::rebind f_;
        BOOST_MPL_PP_REPEAT(i, AUX_LAMBDA_INVOCATION, T)

        typedef mpl::protect< BOOST_PP_CAT(bind,i)<
              f_
            , AUX_LAMBDA_PARAMS(i, arg)
            > > type;
    };
};

#   undef AUX_LAMBDA_INVOCATION

#undef i
#endif // BOOST_PP_IS_ITERATING

//-----------------------------------------------------------------------------
// boost mpl/arg.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Peter Dimov, Aleksey Gurtovoy
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

#ifndef BOOST_MPL_ARG_HPP_INCLUDED
#define BOOST_MPL_ARG_HPP_INCLUDED

#include "boost/mpl/aux_/config/static_constant.hpp"

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/arg_fwd.hpp"
#   include "boost/mpl/void.hpp"
#   include "boost/mpl/aux_/arity_spec.hpp"
#   include "boost/mpl/aux_/arg_typedef.hpp"
#   include "boost/static_assert.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
 && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER arg.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else
#   include "boost/mpl/limits/arity.hpp"
#   include "boost/mpl/aux_/preprocessor/default_params.hpp"
#   include "boost/mpl/aux_/preprocessor/params.hpp"
#   include "boost/mpl/aux_/config/lambda.hpp"
#   include "boost/mpl/aux_/config/dtp.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"

#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

namespace boost {
namespace mpl {

// local macro, #undef-ined at the end of the header
#if !defined(BOOST_NO_DEFAULT_TEMPLATE_PARAMETERS_IN_NESTED_TEMPLATES)
#   define AUX_ARG_N_DEFAULT_PARAMS(param,value) \
    BOOST_MPL_PP_DEFAULT_PARAMS( \
          BOOST_MPL_METAFUNCTION_MAX_ARITY \
        , param \
        , value \
        ) \
    /**/
#else
#   define AUX_ARG_N_DEFAULT_PARAMS(param,value) \
    BOOST_MPL_PP_PARAMS( \
          BOOST_MPL_METAFUNCTION_MAX_ARITY \
        , param \
        ) \
    /**/
#endif

#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_METAFUNCTION_MAX_ARITY, "boost/mpl/arg.hpp"))
#include BOOST_PP_ITERATE()


#   undef AUX_ARG_N_DEFAULT_PARAMS

BOOST_MPL_AUX_NONTYPE_ARITY_SPEC(1,int,arg)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_ARG_HPP_INCLUDED

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

#if i > 0

template<> struct arg<i>
{
    BOOST_STATIC_CONSTANT(int, value = i);
    typedef arg<BOOST_PP_INC(i)> next;
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          AUX_ARG_N_DEFAULT_PARAMS(typename U, void_)
        >
    struct apply
    {
        typedef BOOST_PP_CAT(U,i) type;
#if !defined(__BORLANDC__) || (__BORLANDC__ > 0x561 && defined(BOOST_STRICT_CONFIG))
     private:
        BOOST_STATIC_CONSTANT(bool, nv = !is_void_<type>::value);
        BOOST_STATIC_ASSERT(nv);
#endif
    };
};

#else

template<> struct arg<-1>
{
    BOOST_STATIC_CONSTANT(int, value = -1);
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          AUX_ARG_N_DEFAULT_PARAMS(typename U, void_)
        >
    struct apply
    {
        typedef U1 type;
#if !defined(__BORLANDC__) || (__BORLANDC__ > 0x561 && defined(BOOST_STRICT_CONFIG))
     private:
        BOOST_STATIC_CONSTANT(bool, nv = !is_void_<type>::value);
        BOOST_STATIC_ASSERT(nv);
#endif
    };
};

#endif // i > 0

#undef i
#endif // BOOST_PP_IS_ITERATING

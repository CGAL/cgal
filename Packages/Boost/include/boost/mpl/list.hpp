
#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_LIST_HPP_INCLUDED
#define BOOST_MPL_LIST_HPP_INCLUDED

// Copyright (c) 2000-04 Aleksey Gurtovoy
//
// Use, modification and distribution are subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy 
// at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Source$
// $Date$
// $Revision$

// no include guards, the header is intended for multiple inclusion!

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/limits/list.hpp"
#   include "boost/mpl/void.hpp"
#   include "boost/mpl/aux_/config/preprocessor.hpp"

#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"
#   include "boost/preprocessor/stringize.hpp"

#if !defined(BOOST_NEEDS_TOKEN_PASTING_OP_FOR_TOKENS_JUXTAPOSING)
#   define MPL_AUX_LIST_HEADER \
    BOOST_PP_CAT(list,BOOST_MPL_LIMIT_LIST_SIZE).hpp \
    /**/
#else
#   define MPL_AUX_LIST_HEADER \
    BOOST_PP_CAT(list,BOOST_MPL_LIMIT_LIST_SIZE)##.hpp \
    /**/
#endif

#   include BOOST_PP_STRINGIZE(boost/mpl/list/MPL_AUX_LIST_HEADER)
#   undef MPL_AUX_LIST_HEADER
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if    !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER list.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/aux_/config/nttp.hpp"
#   include "boost/mpl/limits/list.hpp"

#   include "boost/preprocessor/arithmetic/sub.hpp"
#   include "boost/preprocessor/tuple/elem.hpp"
#   include "boost/preprocessor/enum_params_with_a_default.hpp"
#   include "boost/preprocessor/enum_params.hpp"
#   include "boost/preprocessor/enum.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/expr_if.hpp"
#   include "boost/preprocessor/iterate.hpp"

#   include "boost/config.hpp"

namespace boost {
namespace mpl {

#   define AUX_LIST_PARAMS(param) \
    BOOST_PP_ENUM_PARAMS( \
          BOOST_MPL_LIMIT_LIST_SIZE \
        , param \
        ) \
    /**/

#   define AUX_LIST_DEFAULT_PARAMS(param, value) \
     BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT( \
          BOOST_MPL_LIMIT_LIST_SIZE \
        , param \
        , value \
        ) \
    /**/

#   define AUX_LIST_N_PARAMS(n, param) \
    BOOST_PP_ENUM_PARAMS(n, param) \
    /**/

#   define AUX_LIST_N_PARTIAL_SPEC_PARAMS(n, param, def) \
    BOOST_PP_ENUM_PARAMS(n, param) \
    BOOST_PP_COMMA_IF(n) \
    BOOST_PP_ENUM( \
          BOOST_PP_SUB_D(1,BOOST_MPL_LIMIT_LIST_SIZE,n) \
        , BOOST_PP_TUPLE_ELEM_3_2 \
        , def \
        ) \
    /**/

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
// forward declaration
template<
      AUX_LIST_DEFAULT_PARAMS(typename T, void_)
    >
struct list;
#else
namespace aux {
template< BOOST_MPL_AUX_NTTP_DECL(int, N) > struct list_impl_chooser;
}
#endif

#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_LIMIT_LIST_SIZE, "boost/mpl/list.hpp"))
#include BOOST_PP_ITERATE()

// real C++ version is already taken care of
#if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace aux {
// list_count_args
#define BOOST_MPL_AUX_COUNT_ARGS_PREFIX list
#define BOOST_MPL_AUX_COUNT_ARGS_DEFAULT void_
#define BOOST_MPL_AUX_COUNT_ARGS_ARITY BOOST_MPL_LIMIT_LIST_SIZE
#define BOOST_MPL_AUX_COUNT_ARGS_USE_STANDARD_PP_PRIMITIVES
#include "boost/mpl/aux_/count_args.hpp"

template<
      AUX_LIST_PARAMS(typename T)
    >
struct list_impl
{
    typedef aux::list_count_args< AUX_LIST_PARAMS(T) > arg_num_;
    typedef typename aux::list_impl_chooser< 
          arg_num_::value
        >::template result_< AUX_LIST_PARAMS(T) >::type type;
};

} // namespace aux

template<
      AUX_LIST_DEFAULT_PARAMS(typename T, void_)
    >
struct list
    : aux::list_impl< AUX_LIST_PARAMS(T) >::type
{
    typedef typename aux::list_impl<
          AUX_LIST_PARAMS(T)
        >::type type;
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#   undef AUX_LIST_N_PARTIAL_SPEC_PARAMS
#   undef AUX_LIST_N_PARAMS
#   undef AUX_LIST_DEFAULT_PARAMS
#   undef AUX_LIST_PARAMS

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_LIST_HPP_INCLUDED

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

#   if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#if i == BOOST_MPL_LIMIT_LIST_SIZE

//: primary template (not a specialization!)
template<
      AUX_LIST_N_PARAMS(i, typename T)
    >
struct list
    : BOOST_PP_CAT(list,i)< AUX_LIST_N_PARAMS(i, T) >
{
    typedef BOOST_PP_CAT(list,i)< AUX_LIST_N_PARAMS(i, T) > type;
};

#else

template<
      AUX_LIST_N_PARAMS(i, typename T)
    >
struct list< AUX_LIST_N_PARTIAL_SPEC_PARAMS(i, T, void_) >
    : BOOST_PP_CAT(list,i)< AUX_LIST_N_PARAMS(i, T) >
{
    typedef BOOST_PP_EXPR_IF(i,typename) BOOST_PP_CAT(list,i)< AUX_LIST_N_PARAMS(i, T) >::type type;
};

#endif // i == BOOST_MPL_LIMIT_LIST_SIZE

#   else

namespace aux {

template<>
struct list_impl_chooser<i>
{
    template<
          AUX_LIST_PARAMS(typename T)
        >
    struct result_
    {
        typedef typename BOOST_PP_CAT(list,i)<
              AUX_LIST_N_PARAMS(i, T)
            >::type type;
    };
};

} // namespace aux

#   endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#undef i
#endif // BOOST_PP_IS_ITERATING

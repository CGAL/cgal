
#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_LIST_C_HPP_INCLUDED
#define BOOST_MPL_LIST_C_HPP_INCLUDED

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
#   include "boost/mpl/aux_/config/preprocessor.hpp"

#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"
#   include "boost/preprocessor/stringize.hpp"

#if !defined(BOOST_NEEDS_TOKEN_PASTING_OP_FOR_TOKENS_JUXTAPOSING)
#   define MPL_AUX_LIST_HEADER \
    BOOST_PP_CAT(BOOST_PP_CAT(list,BOOST_MPL_LIMIT_LIST_SIZE), _c).hpp \
    /**/
#else
#   define MPL_AUX_LIST_HEADER \
    BOOST_PP_CAT(BOOST_PP_CAT(list,BOOST_MPL_LIMIT_LIST_SIZE), _c)##.hpp \
    /**/
#endif

#   include BOOST_PP_STRINGIZE(boost/mpl/list/MPL_AUX_LIST_HEADER)
#   undef MPL_AUX_LIST_HEADER
#   include <climits>
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) && \
    !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER list_c.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/limits/list.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"

#   include "boost/preprocessor/arithmetic/sub.hpp"
#   include "boost/preprocessor/tuple/elem.hpp"
#   include "boost/preprocessor/enum_params_with_a_default.hpp"
#   include "boost/preprocessor/enum_params.hpp"
#   include "boost/preprocessor/enum.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/iterate.hpp"

#   include "boost/config.hpp"

#if defined(BOOST_MPL_PREPROCESSING_MODE)
#   undef LONG_MAX
#endif

namespace boost {
namespace mpl {

#   define AUX_LIST_C(i) \
    BOOST_PP_CAT(BOOST_PP_CAT(list,i),_c) \
    /**/

#   define AUX_LIST_C_PARAMS(param) \
    BOOST_PP_ENUM_PARAMS( \
          BOOST_MPL_LIMIT_LIST_SIZE \
        , param \
        ) \
    /**/

#   define AUX_LIST_C_DEFAULT_PARAMS(param, value) \
     BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT( \
          BOOST_MPL_LIMIT_LIST_SIZE \
        , param \
        , value \
        ) \
    /**/

#   define AUX_LIST_C_N_PARAMS(n, param) \
    BOOST_PP_COMMA_IF(n) \
    BOOST_PP_ENUM_PARAMS(n, param) \
    /**/

#   define AUX_LIST_C_N_PARTIAL_SPEC_PARAMS(n, param, def) \
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
      typename T
    , AUX_LIST_C_DEFAULT_PARAMS(BOOST_MPL_AUX_NTTP_DECL(long, C), LONG_MAX)
    >
struct list_c;
#else
namespace aux {
template< BOOST_MPL_AUX_NTTP_DECL(int, N) > struct list_c_impl_chooser;
}
#endif

#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_LIMIT_LIST_SIZE, "boost/mpl/list_c.hpp"))
#include BOOST_PP_ITERATE()

// real C++ version is already taken care of
#if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace aux {
// list_count_args
#define BOOST_MPL_AUX_COUNT_ARGS_PREFIX list_c
#define BOOST_MPL_AUX_COUNT_ARGS_DEFAULT LONG_MAX
#define BOOST_MPL_AUX_COUNT_ARGS_TEMPLATE_PARAM BOOST_MPL_AUX_NTTP_DECL(long, T)
#define BOOST_MPL_AUX_COUNT_ARGS_ARITY BOOST_MPL_LIMIT_LIST_SIZE
#define BOOST_MPL_AUX_COUNT_ARGS_USE_STANDARD_PP_PRIMITIVES
#include "boost/mpl/aux_/count_args.hpp"

template<
      typename T
    , AUX_LIST_C_PARAMS(BOOST_MPL_AUX_NTTP_DECL(long, C))
    >
struct list_c_impl
{
    typedef aux::list_c_count_args< AUX_LIST_C_PARAMS(C) > arg_num_;
    typedef typename aux::list_c_impl_chooser< arg_num_::value >
        ::template result_< T, AUX_LIST_C_PARAMS(C) >::type type;
};

} // namespace aux

template<
      typename T
    , AUX_LIST_C_DEFAULT_PARAMS(BOOST_MPL_AUX_NTTP_DECL(long, C), LONG_MAX)
    >
struct list_c
    : aux::list_c_impl< T,AUX_LIST_C_PARAMS(C) >::type
{
    typedef typename aux::list_c_impl<
          T,AUX_LIST_C_PARAMS(C)
        >::type type;
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#   undef AUX_LIST_C_N_PARTIAL_SPEC_PARAMS
#   undef AUX_LIST_C_N_PARAMS
#   undef AUX_LIST_C_DEFAULT_PARAMS
#   undef AUX_LIST_C_PARAMS
#   undef AUX_LIST_C

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_LIST_C_HPP_INCLUDED

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

#   if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#if i == BOOST_MPL_LIMIT_LIST_SIZE

//: primary template (not a specialization!)
template<
      typename T
    AUX_LIST_C_N_PARAMS(i, BOOST_MPL_AUX_NTTP_DECL(long, C))
    >
struct list_c
    : AUX_LIST_C(i)< T AUX_LIST_C_N_PARAMS(i, C) >
{
    typedef typename AUX_LIST_C(i)< T AUX_LIST_C_N_PARAMS(i, C) >::type type;
};

#else

template<
      typename T
    AUX_LIST_C_N_PARAMS(i, BOOST_MPL_AUX_NTTP_DECL(long, C))
    >
struct list_c< T,AUX_LIST_C_N_PARTIAL_SPEC_PARAMS(i, C, LONG_MAX) >
    : AUX_LIST_C(i)< T AUX_LIST_C_N_PARAMS(i, C) >
{
    typedef typename AUX_LIST_C(i)< T AUX_LIST_C_N_PARAMS(i, C) >::type type;
};

#endif // i == BOOST_MPL_LIMIT_LIST_SIZE

#   else

namespace aux {

template<>
struct list_c_impl_chooser<i>
{
    template<
          typename T
        , AUX_LIST_C_PARAMS(BOOST_MPL_AUX_NTTP_DECL(long, C))
        >
    struct result_
    {
        typedef typename AUX_LIST_C(i)<
              T AUX_LIST_C_N_PARAMS(i, C)
            >::type type;
    };
};

} // namespace aux

#   endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#undef i
#endif // BOOST_PP_IS_ITERATING

//-----------------------------------------------------------------------------
// boost mpl/aux_/count_args.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

// no include guards, the header is intended for multiple inclusion!

#include "boost/preprocessor/expr_if.hpp"
#include "boost/preprocessor/inc.hpp"
#include "boost/preprocessor/cat.hpp"

#if !defined(BOOST_MPL_AUX_COUNT_ARGS_TEMPLATE_PARAM)
#   define BOOST_MPL_AUX_COUNT_ARGS_TEMPLATE_PARAM typename T
#endif

// local macros, #undef-ined at the end of the header

#if !defined(BOOST_MPL_AUX_COUNT_ARGS_USE_STANDARD_PP_PRIMITIVES)

#   include "boost/mpl/aux_/preprocessor/repeat.hpp"
#   include "boost/mpl/aux_/preprocessor/params.hpp"

#   define AUX_COUNT_ARGS_REPEAT BOOST_MPL_PP_REPEAT
#   define AUX_COUNT_ARGS_PARAMS(param) \
    BOOST_MPL_PP_PARAMS( \
          BOOST_MPL_AUX_COUNT_ARGS_ARITY \
        , param \
        ) \
    /**/

#else

#   include "boost/preprocessor/enum_shifted_params.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/inc.hpp"

#   define AUX_COUNT_ARGS_REPEAT BOOST_PP_REPEAT_1
#   define AUX_COUNT_ARGS_PARAMS(param) \
    BOOST_PP_ENUM_SHIFTED_PARAMS( \
          BOOST_PP_INC(BOOST_MPL_AUX_COUNT_ARGS_ARITY) \
        , param \
        ) \
    /**/

#endif // BOOST_MPL_AUX_COUNT_ARGS_USE_STANDARD_PP_PRIMITIVES


#define AUX_IS_ARG_TEMPLATE_NAME \
    BOOST_PP_CAT(is_,BOOST_PP_CAT(BOOST_MPL_AUX_COUNT_ARGS_PREFIX,_arg)) \
/**/

#define AUX_COUNT_ARGS_FUNC(unused, i, param) \
    BOOST_PP_EXPR_IF(i, +) \
    AUX_IS_ARG_TEMPLATE_NAME<BOOST_PP_CAT(param,BOOST_PP_INC(i))>::value \
/**/

// is_<xxx>_arg
template< BOOST_MPL_AUX_COUNT_ARGS_TEMPLATE_PARAM >
struct AUX_IS_ARG_TEMPLATE_NAME
{
    BOOST_STATIC_CONSTANT(bool, value = true);
};

template<>
struct AUX_IS_ARG_TEMPLATE_NAME<BOOST_MPL_AUX_COUNT_ARGS_DEFAULT>
{
    BOOST_STATIC_CONSTANT(bool, value = false);
};

// <xxx>_count_args
template<
      AUX_COUNT_ARGS_PARAMS(BOOST_MPL_AUX_COUNT_ARGS_TEMPLATE_PARAM)
    >
struct BOOST_PP_CAT(BOOST_MPL_AUX_COUNT_ARGS_PREFIX,_count_args)
{
    BOOST_STATIC_CONSTANT(int, value = AUX_COUNT_ARGS_REPEAT(
          BOOST_MPL_AUX_COUNT_ARGS_ARITY
        , AUX_COUNT_ARGS_FUNC
        , T
        ));
};

#undef AUX_COUNT_ARGS_FUNC
#undef AUX_IS_ARG_TEMPLATE_NAME
#undef AUX_COUNT_ARGS_PARAMS
#undef AUX_COUNT_ARGS_REPEAT

#undef BOOST_MPL_AUX_COUNT_ARGS_ARITY
#undef BOOST_MPL_AUX_COUNT_ARGS_DEFAULT
#undef BOOST_MPL_AUX_COUNT_ARGS_PREFIX
#undef BOOST_MPL_AUX_COUNT_ARGS_USE_STANDARD_PP_PRIMITIVES
#undef BOOST_MPL_AUX_COUNT_ARGS_TEMPLATE_PARAM

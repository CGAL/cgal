
// + file: boost/mpl/aux_/logical_op.hpp
// + last modified: 25/feb/03

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
//
// See http://www.boost.org/libs/mpl for documentation.

// no include guards, the header is intended for multiple inclusion!

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/bool.hpp"
#   include "boost/mpl/aux_/nested_type_wknd.hpp"
#   include "boost/mpl/aux_/void_spec.hpp"
#   include "boost/mpl/aux_/lambda_support.hpp"
#endif

#include "boost/mpl/limits/arity.hpp"
#include "boost/mpl/aux_/preprocessor/params.hpp"
#include "boost/mpl/aux_/preprocessor/ext_params.hpp"
#include "boost/mpl/aux_/preprocessor/def_params_tail.hpp"
#include "boost/mpl/aux_/preprocessor/enum.hpp"
#include "boost/mpl/aux_/preprocessor/sub.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

#include "boost/preprocessor/dec.hpp"
#include "boost/preprocessor/inc.hpp"
#include "boost/preprocessor/cat.hpp"

namespace boost { namespace mpl {

// local macros, #undef-ined at the end of the header
#   define AUX_LOGICAL_OP_PARAMS(param, sub) \
    BOOST_MPL_PP_PARAMS( \
          BOOST_MPL_PP_SUB(BOOST_MPL_METAFUNCTION_MAX_ARITY, sub) \
        , param \
        ) \
    /**/

#   define AUX_LOGICAL_OP_SHIFTED_PARAMS(param, sub) \
    BOOST_MPL_PP_EXT_PARAMS( \
          2, BOOST_MPL_PP_SUB(BOOST_PP_INC(BOOST_MPL_METAFUNCTION_MAX_ARITY), sub) \
        , param \
        ) \
    /**/

#   define AUX_LOGICAL_OP_SPEC_PARAMS(param) \
    BOOST_MPL_PP_ENUM( \
          BOOST_PP_DEC(BOOST_MPL_METAFUNCTION_MAX_ARITY) \
        , param \
        ) \
    /**/

namespace aux {

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template< bool C_, AUX_LOGICAL_OP_PARAMS(typename T, 1) >
struct BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)
    : BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE1,_)
{
};

template< AUX_LOGICAL_OP_PARAMS(typename T, 1) >
struct BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)< AUX_LOGICAL_OP_VALUE2,AUX_LOGICAL_OP_PARAMS(T, 1) >
    : BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)<
          BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
        , AUX_LOGICAL_OP_SHIFTED_PARAMS(T, 1)
        , BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_)
        >
{
};

template<>
struct BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)<
          AUX_LOGICAL_OP_VALUE2
        , AUX_LOGICAL_OP_SPEC_PARAMS(BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_))
        >
    : BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_)
{
};

#else

template< bool C_ > struct BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)
{
    template< AUX_LOGICAL_OP_PARAMS(typename T, 1) > struct result_
        : BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE1,_)
    {
    };
};

template<> struct BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)<AUX_LOGICAL_OP_VALUE2>
{
    template< AUX_LOGICAL_OP_PARAMS(typename T, 1) > struct result_
        : BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)< 
              BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
            >::template result_< AUX_LOGICAL_OP_SHIFTED_PARAMS(T,1),BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_) >
    {
    };

#if BOOST_WORKAROUND(BOOST_MSVC, == 1300)
    template<> struct result_<AUX_LOGICAL_OP_SPEC_PARAMS(BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_))>
        : BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_)
    {
    };
};
#else
};

template<>
struct BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)<AUX_LOGICAL_OP_VALUE2>
    ::result_< AUX_LOGICAL_OP_SPEC_PARAMS(BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_)) >
        : BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_)
{
};
#endif // BOOST_MSVC == 1300

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

} // namespace aux

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    BOOST_MPL_PP_DEF_PARAMS_TAIL(2, typename T, BOOST_PP_CAT(AUX_LOGICAL_OP_VALUE2,_))
    >
struct AUX_LOGICAL_OP_NAME
#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
    : aux::BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)<
          BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
        , AUX_LOGICAL_OP_SHIFTED_PARAMS(T,0)
        >
#else
    : aux::BOOST_PP_CAT(AUX_LOGICAL_OP_NAME,impl)< 
          BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
        >::template result_< AUX_LOGICAL_OP_SHIFTED_PARAMS(T,0) >
#endif
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          BOOST_MPL_METAFUNCTION_MAX_ARITY
        , AUX_LOGICAL_OP_NAME
        , (AUX_LOGICAL_OP_PARAMS(T, 0))
        )
};

BOOST_MPL_AUX_VOID_SPEC_EXT(
      2
    , BOOST_MPL_METAFUNCTION_MAX_ARITY
    , AUX_LOGICAL_OP_NAME
    )

}} // namespace boost::mpl

#undef AUX_LOGICAL_OP_SPEC_PARAMS
#undef AUX_LOGICAL_OP_SHIFTED_PARAMS
#undef AUX_LOGICAL_OP_PARAMS
#undef AUX_LOGICAL_OP_NAME
#undef AUX_LOGICAL_OP_VALUE1
#undef AUX_LOGICAL_OP_VALUE2

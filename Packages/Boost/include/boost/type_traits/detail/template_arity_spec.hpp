//-----------------------------------------------------------------------------
// boost/type_traits/detail/template_arity_spec.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Use, modification and distribution are subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

// no include guards, the header is intended for multiple inclusion!

#include "boost/mpl/aux_/template_arity_fwd.hpp"
#include "boost/mpl/aux_/preprocessor/params.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/mpl/aux_/config/overload_resolution.hpp"
#include "boost/config.hpp"

#if defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT) && \
    defined(BOOST_MPL_BROKEN_OVERLOAD_RESOLUTION)
#   define BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(i, name) \
namespace mpl { namespace aux { \
template< BOOST_MPL_PP_PARAMS(i, typename T) > \
struct template_arity< \
      name< BOOST_MPL_PP_PARAMS(i, T) > \
    > \
{ \
    BOOST_STATIC_CONSTANT(int, value = i ); \
}; \
}} \
/**/
#else
#   define BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(i, name) /**/
#endif

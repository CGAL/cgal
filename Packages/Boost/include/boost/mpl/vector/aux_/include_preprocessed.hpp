
// Copyright (c) 2001-04 Aleksey Gurtovoy
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

#include "boost/mpl/aux_/config/vector.hpp"
#include "boost/mpl/aux_/config/ctps.hpp"
#include "boost/mpl/aux_/config/preprocessor.hpp"
#include "boost/preprocessor/cat.hpp"
#include "boost/preprocessor/stringize.hpp"

#if defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL)
#   define AUX_VECTOR_INCLIDE_DIR typeof_based
#elif defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
   || defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION)
#   define AUX_VECTOR_INCLIDE_DIR no_ctps
#else
#   define AUX_VECTOR_INCLIDE_DIR plain
#endif

#if !defined(BOOST_NEEDS_TOKEN_PASTING_OP_FOR_TOKENS_JUXTAPOSING)
#   define AUX_PREPROCESSED_HEADER \
    AUX_VECTOR_INCLIDE_DIR/BOOST_MPL_PREPROCESSED_HEADER \
/**/
#else
#   define AUX_PREPROCESSED_HEADER \
    BOOST_PP_CAT(AUX_VECTOR_INCLIDE_DIR,/)##BOOST_MPL_PREPROCESSED_HEADER \
/**/
#endif


#   include BOOST_PP_STRINGIZE(boost/mpl/vector/aux_/preprocessed/AUX_PREPROCESSED_HEADER)

#   undef AUX_PREPROCESSED_HEADER
#   undef AUX_VECTOR_INCLIDE_DIR

#undef BOOST_MPL_PREPROCESSED_HEADER

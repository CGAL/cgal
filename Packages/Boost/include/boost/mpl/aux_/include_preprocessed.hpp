
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

#include "boost/mpl/aux_/config/compiler.hpp"
#include "boost/mpl/aux_/config/preprocessor.hpp"
#include "boost/preprocessor/cat.hpp"
#include "boost/preprocessor/stringize.hpp"

#if !defined(BOOST_NEEDS_TOKEN_PASTING_OP_FOR_TOKENS_JUXTAPOSING)
#   define AUX_PREPROCESSED_HEADER \
    BOOST_MPL_COMPILER_DIR/BOOST_MPL_PREPROCESSED_HEADER \
/**/
#else
#   define AUX_PREPROCESSED_HEADER \
    BOOST_PP_CAT(BOOST_MPL_COMPILER_DIR,/)##BOOST_MPL_PREPROCESSED_HEADER \
/**/
#endif

#   include BOOST_PP_STRINGIZE(boost/mpl/aux_/preprocessed/AUX_PREPROCESSED_HEADER)
#   undef AUX_PREPROCESSED_HEADER

#undef BOOST_MPL_PREPROCESSED_HEADER


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

#include "boost/preprocessor/cat.hpp"
#include "boost/preprocessor/stringize.hpp"

#   define AUX_PREPROCESSED_HEADER \
    aux_/preprocessed/plain/BOOST_MPL_PREPROCESSED_HEADER \
/**/

#   include BOOST_PP_STRINGIZE(boost/mpl/list/AUX_PREPROCESSED_HEADER)

#   undef AUX_PREPROCESSED_HEADER

#undef BOOST_MPL_PREPROCESSED_HEADER


#ifndef BOOST_MPL_AND_HPP_INCLUDED
#define BOOST_MPL_AND_HPP_INCLUDED

// + file: boost/mpl/and.hpp
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

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
 && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   include "boost/mpl/bool.hpp"
#   include "boost/mpl/aux_/nested_type_wknd.hpp"
#   include "boost/mpl/aux_/void_spec.hpp"
#   include "boost/mpl/aux_/lambda_support.hpp"

#   define BOOST_MPL_PREPROCESSED_HEADER and.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   define AUX_LOGICAL_OP_NAME and_
#   define AUX_LOGICAL_OP_VALUE1 false
#   define AUX_LOGICAL_OP_VALUE2 true
#   include "boost/mpl/aux_/logical_op.hpp"

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AND_HPP_INCLUDED


#ifndef BOOST_MPL_SIZE_T_HPP_INCLUDED
#define BOOST_MPL_SIZE_T_HPP_INCLUDED

// + file: boost/mpl/size_t.hpp
// + last modified: 08/mar/03

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

#include "boost/mpl/size_t_fwd.hpp"

#define AUX_WRAPPER_VALUE_TYPE std::size_t
#define AUX_WRAPPER_NAME size_t
#define AUX_WRAPPER_PARAMS(N) std::size_t N

#include "boost/mpl/aux_/integral_wrapper.hpp"

#endif // BOOST_MPL_SIZE_T_HPP_INCLUDED

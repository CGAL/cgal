
#ifndef BOOST_MPL_COUNT_FWD_HPP_INCLUDED
#define BOOST_MPL_COUNT_FWD_HPP_INCLUDED

// + file: boost/mpl/count_fwd.hpp
// + last modified: 05/nov/03

// Copyright Aleksey Gurtovoy 2000-03
//
// Use, modification and distribution are subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy 
// at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/aux_/algorithm_namespace.hpp"

namespace boost {
namespace mpl {

template< typename Tag > struct count_impl;

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN
template< typename Sequence, typename T > struct count;
BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

}}

#endif // BOOST_MPL_COUNT_FWD_HPP_INCLUDED

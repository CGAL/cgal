//-----------------------------------------------------------------------------
// boost mpl/vector/aux_/tag.hpp header file
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

#ifndef BOOST_MPL_VECTOR_AUX_TAG_HPP_INCLUDED
#define BOOST_MPL_VECTOR_AUX_TAG_HPP_INCLUDED

#include "boost/mpl/aux_/config/vector.hpp"

namespace boost {
namespace mpl {
namespace aux {

#if defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL)
struct vector_tag;
#else
template< long N > struct vector_tag;
#endif

} // namespace aux
} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_VECTOR_AUX_TAG_HPP_INCLUDED

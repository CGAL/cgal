
#ifndef BOOST_MPL_INSERT_FWD_HPP_INCLUDED
#define BOOST_MPL_INSERT_FWD_HPP_INCLUDED

// + file: boost/mpl/insert_fwd.hpp
// + last modified: 02/may/03

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

namespace boost {
namespace mpl {

template< typename Tag > struct insert_traits;
template< typename Sequence, typename Pos_or_T, typename T > struct insert;

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_INSERT_FWD_HPP_INCLUDED

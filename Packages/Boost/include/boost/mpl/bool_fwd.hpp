
#ifndef BOOST_MPL_BOOL_FWD_HPP_INCLUDED
#define BOOST_MPL_BOOL_FWD_HPP_INCLUDED

// + file: boost/mpl/bool_fwd.hpp
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

namespace boost { namespace mpl {

template< bool C_ > struct bool_;

// shorcuts
typedef bool_<true> true_;
typedef bool_<false> false_;

}}

#endif // BOOST_MPL_BOOL_FWD_HPP_INCLUDED

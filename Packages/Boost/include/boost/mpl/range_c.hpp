//-----------------------------------------------------------------------------
// boost mpl/range_c.hpp header file
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

#ifndef BOOST_MPL_RANGE_C_HPP_INCLUDED
#define BOOST_MPL_RANGE_C_HPP_INCLUDED

#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/aux_/range_c/front.hpp"
#include "boost/mpl/aux_/range_c/back.hpp"
#include "boost/mpl/aux_/range_c/size.hpp"
#include "boost/mpl/aux_/range_c/O1_size.hpp"
#include "boost/mpl/aux_/range_c/empty.hpp"
#include "boost/mpl/aux_/range_c/iterator.hpp"
#include "boost/mpl/aux_/range_c/tag.hpp"

namespace boost {
namespace mpl {

template<
      typename T
    , T Start
    , T Finish
    >
struct range_c
{
    typedef aux::half_open_range_tag tag;
    typedef range_c type;

    typedef integral_c<T,Start> start;
    typedef integral_c<T,Finish> finish;

    typedef range_c_iterator<start> begin;
    typedef range_c_iterator<finish> end;
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_RANGE_C_HPP_INCLUDED

//-----------------------------------------------------------------------------
// boost mpl/sequence_tag_fwd.hpp header file
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

#ifndef BOOST_MPL_SEQUENCE_TAG_FWD_HPP_INCLUDED
#define BOOST_MPL_SEQUENCE_TAG_FWD_HPP_INCLUDED

namespace boost { namespace mpl {

struct nested_begin_end_tag;
struct non_sequence_tag;

template< typename Sequence > struct sequence_tag;

}} // namespace boost::mpl

#endif // BOOST_MPL_SEQUENCE_TAG_FWD_HPP_INCLUDED

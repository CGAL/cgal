//-----------------------------------------------------------------------------
// boost mpl/iterator_tag.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_ITERATOR_TAG_HPP_INCLUDED
#define BOOST_MPL_ITERATOR_TAG_HPP_INCLUDED

namespace boost {
namespace mpl {

struct input_iter_tag_;
struct fwd_iter_tag_;
struct bi_iter_tag_;
struct ra_iter_tag_;

typedef input_iter_tag_ input_iterator_tag;
typedef fwd_iter_tag_   forward_iterator_tag;
typedef bi_iter_tag_    bidirectional_iterator_tag;
typedef ra_iter_tag_    random_access_iterator_tag;

} // namespace mpl
} // namespace boost 

#endif // BOOST_MPL_ITERATOR_TAG_HPP_INCLUDED

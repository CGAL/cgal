// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL_TAGS_H
#define CGAL_TAGS_H

#include <CGAL/IO/io_tags.h>
#include <boost/mpl/integral_c.hpp>

namespace CGAL {

struct Void {};

// Boolean_tag<bool> is a model of the Boost Integral Constant concept.
// https://www.boost.org/libs/mpl/doc/refmanual/integral-constant.html
template <bool b>
struct Boolean_tag {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  static const bool value = b;
  typedef Boolean_tag<b> type;
  operator bool() const { return this->value; }
};
/* In C++11, try:
template <bool b>
using Boolean_tag = std::integral_constant<bool, b>;
*/

typedef Boolean_tag<true>   Tag_true;
typedef Boolean_tag<false>  Tag_false;

// the function check_tag is deprecated since CGAL 3.3
inline bool check_tag( Tag_true)  {return true;}
inline bool check_tag( Tag_false) {return false;}

struct Null_tag {};

struct Null_functor {
  typedef Null_tag result_type;
  typedef Null_tag second_argument_type;
};

// For concurrency
struct Sequential_tag {};
struct Parallel_tag : public Sequential_tag {};

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Parallel_if_available_tag;
#else
typedef CGAL::Sequential_tag Parallel_if_available_tag;
#endif

// A function that asserts a specific compile time tag
// forcing its two arguments to have equal type.
template <class Base>
struct Assert_tag_class
{
    void match_compile_time_tag( const Base&) const {}
};

template <class Tag, class Derived>
inline
void
Assert_compile_time_tag( const Tag&, const Derived& b)
{
  Assert_tag_class<Tag> x;
  x.match_compile_time_tag(b);
}

} //namespace CGAL

#endif // CGAL_TAGS_H

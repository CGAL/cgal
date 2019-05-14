// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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
//   http://www.boost.org/libs/mpl/doc/refmanual/integral-constant.html
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

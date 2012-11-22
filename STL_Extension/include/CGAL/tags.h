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
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_TAGS_H
#define CGAL_TAGS_H

#include <CGAL/IO/io_tags.h>
#include <boost/config.hpp>

namespace CGAL {

struct Void {};

// Two struct's to denote boolean compile time decisions.

template <bool b>
struct Boolean_tag {
  static const bool value = b;
};

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

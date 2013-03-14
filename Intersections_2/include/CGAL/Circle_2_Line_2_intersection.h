// Copyright (c) 2000  
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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_CIRCLE_2_LINE_2_INTERSECTION_H
#define CGAL_CIRCLE_2_LINE_2_INTERSECTION_H

#include <CGAL/Circle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/squared_distance_2_1.h>

namespace CGAL {

namespace internal {

template <class K>
bool
do_intersect(const typename K::Circle_2 & c, 
	     const typename K::Line_2& l,
	     const K&)
{
    return squared_distance(c.center(), l) <= c.squared_radius();
}

template <class K>
bool
do_intersect(const typename K::Line_2& l, 
	     const typename K::Circle_2 & c,
	     const K&)
{
    return squared_distance(c.center(), l) <= c.squared_radius();
}

} // namespace internal

template <class K>
inline
bool
do_intersect(const Circle_2<K> & c, 
	     const Line_2<K> & l)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(c, l);
}

template <class K>
inline
bool
do_intersect(const Line_2<K> & l, 
	     const Circle_2<K> & c)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(c, l);
}


} //namespace CGAL

#endif

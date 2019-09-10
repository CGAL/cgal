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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_CIRCLE_2_CIRCLE_2_INTERSECTION_H
#define CGAL_CIRCLE_2_CIRCLE_2_INTERSECTION_H

#include <CGAL/Circle_2.h>
#include <CGAL/squared_distance_2_1.h>

namespace CGAL {

namespace internal {

template <class K>
bool
do_intersect(const typename K::Circle_2 & circ1, 
	     const typename K::Circle_2& circ2,
	     const K&)
{
    typedef typename K::FT FT;
    FT sr1 = circ1.squared_radius();
    FT sr2 = circ2.squared_radius();
    FT squared_dist = squared_distance(circ1.center(), circ2.center());
    FT temp = sr1+sr2-squared_dist;
    return !(FT(4)*sr1*sr2 < temp*temp);
}

} // namespace internal

template <class K>
inline
bool
do_intersect(const Circle_2<K> & circ1, 
	     const Circle_2<K> & circ2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(circ1, circ2);
}


} //namespace CGAL

#endif

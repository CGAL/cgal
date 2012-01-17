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


#ifndef CGAL_TRIANGLE_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_TRIANGLE_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Object.h>
#include <CGAL/Triangle_2_Triangle_2_do_intersect.h>

namespace CGAL {

template <class K>
inline
Object
intersection(const Triangle_2<K> &tr1, 
	     const Triangle_2<K>& tr2)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(tr1, tr2);
}

} //namespace CGAL

#include <CGAL/Intersections_2/Triangle_2_Triangle_2_intersection_impl.h>

#endif

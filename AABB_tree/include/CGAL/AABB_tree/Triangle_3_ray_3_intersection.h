// Copyrigth (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     :  Laurent Rineau

#ifndef CGAL_TRIANGLE_3_RAY_3_INTERSECTION_H
#define CGAL_TRIANGLE_3_RAY_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>

namespace CGAL {
namespace CGALi {
template <class K>
Object
intersection(const typename K::Triangle_3  &t, 
	     const typename K::Ray_3 &r, 
	     const K& k)
{
  return CGAL::intersection(t.supporting_plane(), 
                            r);
}

} // end namespace CGALi

template <class K>
inline
Object 
intersection(const Triangle_3<K> &t, const Ray_3<K> &r)
{
  return typename K::Intersect_3()(t, r);
}

} // end namespace CGAL

#endif // CGAL_TRIANGLE_3_RAY_3_INTERSECTION_H

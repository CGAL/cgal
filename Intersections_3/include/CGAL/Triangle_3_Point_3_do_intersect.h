// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s)     : Philippe Guigue

#ifndef CGAL_TRIANGLE_3_POINT_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_POINT_3_DO_INTERSECT_H

namespace CGAL {

namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3 &t, 
		  const typename K::Point_3    &p,
		  const K & k )
{

  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t));

  typedef typename K::Point_3 Point_3;
  
  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  
  typename K::Orientation_3 orientation = 
    k.orientation_3_object();

  typename K::Coplanar_orientation_3 coplanar_orientation = 
    k.coplanar_orientation_3_object();



  const Point_3 & a = vertex_on(t,0);
  const Point_3 & b = vertex_on(t,1);
  const Point_3 & c = vertex_on(t,2);
  

  if (orientation(a,b,c,p) != COPLANAR)
    return false;
  

  const Orientation abp = coplanar_orientation(a,b,p);
  const Orientation bcp = coplanar_orientation(b,c,p);
  
  
  switch ( abp ) {
  case POSITIVE: return  bcp != NEGATIVE 
		   &&   coplanar_orientation(c,a,p) != NEGATIVE ;
  case NEGATIVE: return  bcp != POSITIVE  
		   &&   coplanar_orientation(c,a,p) != POSITIVE ;
  case COLLINEAR:
    switch ( bcp ) {
    case POSITIVE: return  coplanar_orientation(c,a,p) != NEGATIVE ;
    case NEGATIVE: return  coplanar_orientation(c,a,p) != POSITIVE ;
    case COLLINEAR: return true;
    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  default: // should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
  
}


template <class K>
bool do_intersect(const typename K::Point_3    &p,
		  const typename K::Triangle_3 &t, 
		  const K & k )
{
  return do_intersect(t, p, k);
}


} // namespace internal

 
template <class K>
inline bool do_intersect(const Point_3<K> &p, 
			 const Triangle_3<K> &t)
{ 
  return typename K::Do_intersect_3()(t,p);
} 

template <class K>
inline bool do_intersect(const Triangle_3<K> &t, 
			 const Point_3<K> &p)
{ 
  return typename K::Do_intersect_3()(t,p);
} 

/*
template <class K>
inline bool do_intersect(const Point_3<K> &p, 
			 const Triangle_3<K> &t,
			 const K & k)
{ 
  return do_intersect(t,p,k);
} 
*/

} //namespace CGAL

#endif // CGAL_TRIANGLE_3_POINT_3_DO_INTERSECT_H

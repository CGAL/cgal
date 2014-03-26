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

#ifndef CGAL_TRIANGLE_3_LINE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_LINE_3_DO_INTERSECT_H

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

  template <class K>
  class Triangle_3;

  template <class K>
  class Line_3;

namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3 &t, 
		  const typename K::Line_3     &l,
		  const K & k )
{
  
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(l) ) ;
  
  typedef typename K::Point_3 Point_3;
  

  typename K::Construct_point_on_3 point_on = 
    k.construct_point_on_3_object();
  
  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  
  typename K::Orientation_3 orientation = 
    k.orientation_3_object();

  typename K::Coplanar_orientation_3 coplanar_orientation = 
    k.coplanar_orientation_3_object();



  const Point_3 & a = vertex_on(t,0);
  const Point_3 & b = vertex_on(t,1);
  const Point_3 & c = vertex_on(t,2);
  const Point_3 & p = point_on(l,0);
  const Point_3 & q = point_on(l,1);


  
  if ( ( orientation(a,b,c,p) != COPLANAR ) || 
       ( orientation(a,b,c,q) != COPLANAR ) )
    {

      const Orientation pqab = orientation(p,q,a,b);
      const Orientation pqbc = orientation(p,q,b,c);

      
      switch ( pqab ) {
      case POSITIVE: return  pqbc != NEGATIVE  &&  
		       orientation(p,q,c,a) != NEGATIVE ;
      case NEGATIVE: return  pqbc != POSITIVE  &&  
		       orientation(p,q,c,a) != POSITIVE ;
      case COPLANAR:
	switch ( pqbc ) {
	case POSITIVE: return  orientation(p,q,c,a) != NEGATIVE ;
	case NEGATIVE: return  orientation(p,q,c,a) != POSITIVE ;
	case COPLANAR: return true;
	default: // should not happen.
	  CGAL_kernel_assertion(false);
	  return false;
	}
      default: // should not happen.
	CGAL_kernel_assertion(false);
	return false;
      }
    }
  
  // Coplanar case
    
  const Orientation pqa = coplanar_orientation(p,q,a);
  
  return  coplanar_orientation(p,q,b) != pqa 
    ||  coplanar_orientation(p,q,c) != pqa ;  
}


template <class K>
inline
bool do_intersect(const typename K::Line_3     &l,
		  const typename K::Triangle_3 &t, 
		  const K & k )
{
  return do_intersect(t, l, k);
}

} // namespace internal

CGAL_DO_INTERSECT_FUNCTION(Triangle_3, Line_3, 3)

} //namespace CGAL

#endif //CGAL_TRIANGLE_3_LINE_3_DO_INTERSECT_H

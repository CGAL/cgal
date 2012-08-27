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

#ifndef CGAL_TRIANGLE_2_TRIANGLE_2_DO_INTERSECT_H
#define CGAL_TRIANGLE_2_TRIANGLE_2_DO_INTERSECT_H

#include <CGAL/Triangle_2.h>

namespace CGAL {



namespace internal {

template <class K>
bool intersection_test_vertex(const typename K::Point_2 *  P1, 
			      const typename K::Point_2 *  Q1, 
			      const typename K::Point_2 *  R1,
			      const typename K::Point_2 *  P2, 
			      const typename K::Point_2 *  Q2, 
			      const typename K::Point_2 *  R2,
			      const K & k ){
  
  
  CGAL_kernel_precondition( k.orientation_2_object() (*P1,*Q1,*R1)
			    == POSITIVE);
  CGAL_kernel_precondition( k.orientation_2_object() (*P2,*Q2,*R2)
			    == POSITIVE);

  typename K::Orientation_2 orientation = 
    k.orientation_2_object();
  
  if (orientation(*R2,*P2,*Q1) != NEGATIVE) {
    if (orientation(*R2,*Q2,*Q1) != POSITIVE) {
      if (orientation(*P1,*P2,*Q1) == POSITIVE)
	return orientation(*P1,*Q2,*Q1) != POSITIVE;
      return   orientation(*P1,*P2,*R1) != NEGATIVE
	&&  orientation(*Q1,*R1,*P2) != NEGATIVE  ;
    }
    return  orientation(*P1,*Q2,*Q1) != POSITIVE
      && orientation(*R2,*Q2,*R1) != POSITIVE
      &&  orientation(*Q1,*R1,*Q2) != NEGATIVE ;
  }
  
  if (orientation(*R2,*P2,*R1) != NEGATIVE) { 
    if (orientation(*Q1,*R1,*R2) != NEGATIVE)
      return orientation(*P1,*P2,*R1) != NEGATIVE;
    return orientation(*Q1,*R1,*Q2) != NEGATIVE
      &&  orientation(*R2,*R1,*Q2) != NEGATIVE;
  }
  
  return false; 
  
}


template <class K>
bool intersection_test_edge(const typename K::Point_2 *  P1, 
			    const typename K::Point_2 *  Q1, 
			    const typename K::Point_2 *  R1,
			    const typename K::Point_2 *  P2, 
			    const typename K::Point_2 *
                              CGAL_kernel_precondition_code(Q2),
			    const typename K::Point_2 *  R2,
			    const K & k ){
  
  
  CGAL_kernel_precondition( k.orientation_2_object() (*P1,*Q1,*R1)
			    == POSITIVE);
  CGAL_kernel_precondition( k.orientation_2_object() (*P2,*Q2,*R2)
			    == POSITIVE);
  
  typename K::Orientation_2 orientation = 
    k.orientation_2_object();
  
  if (orientation(*R2,*P2,*Q1) != NEGATIVE) { 
    if (orientation(*P1,*P2,*Q1) != NEGATIVE) 
      return orientation(*P1,*Q1,*R2) != NEGATIVE;
    return orientation(*Q1,*R1,*P2) != NEGATIVE
      && orientation(*R1,*P1,*P2) != NEGATIVE;
  } 

  if (orientation(*R2,*P2,*R1) != NEGATIVE) 
    return  orientation(*P1,*P2,*R1) != NEGATIVE 
      && ( orientation(*P1,*R1,*R2) != NEGATIVE
	   || orientation(*Q1,*R1,*R2) != NEGATIVE ) ;
  
  return false; 
  
}


template <class K>
bool do_intersect(const typename K::Triangle_2 &t1, 
		  const typename K::Triangle_2 &t2,
		  const K & k ){
  
  CGAL_kernel_precondition( ! k.is_degenerate_2_object() (t1) );
  CGAL_kernel_precondition( ! k.is_degenerate_2_object() (t2) );
  
  typename K::Construct_vertex_2 vertex_on =
    k.construct_vertex_2_object();

  typename K::Orientation_2 orientation = 
    k.orientation_2_object();


  typedef typename K::Point_2 Point_2;
  
  const Point_2 & P1 = vertex_on(t1,0);
  const Point_2 & Q1 = vertex_on(t1,1);
  const Point_2 & R1 = vertex_on(t1,2);
  const Point_2 & P2 = vertex_on(t2,0);
  const Point_2 & Q2 = vertex_on(t2,1);
  const Point_2 & R2 = vertex_on(t2,2);

  const Point_2 * p1 = &P1;
  const Point_2 * q1 = &Q1;
  const Point_2 * r1 = &R1;

  const Point_2 * p2 = &P2;
  const Point_2 * q2 = &Q2;
  const Point_2 * r2 = &R2;

  if ( orientation(P1,Q1,R1) != POSITIVE ) {
    q1 = &R1;;
    r1 = &Q1;
  }
  
  if ( orientation(P2,Q2,R2) != POSITIVE ) {
    q2 = &R2;
    r2 = &Q2;
  }

  if ( orientation(*p2,*q2,*p1) != NEGATIVE ) {
    if ( orientation(*q2,*r2,*p1) != NEGATIVE ) { 
      if ( orientation(*r2,*p2,*p1) != NEGATIVE ) return true;
      return internal::intersection_test_edge(p1,q1,r1,p2,q2,r2,k);
    } 
    if ( orientation(*r2,*p2,*p1) != NEGATIVE ) 
      return internal::intersection_test_edge(p1,q1,r1,r2,p2,q2,k);
    return internal::intersection_test_vertex(p1,q1,r1,p2,q2,r2,k);
    
  }
  
  if ( orientation(*q2,*r2,*p1) != NEGATIVE ) {
    if ( orientation(*r2,*p2,*p1) != NEGATIVE ) 
      return internal::intersection_test_edge(p1,q1,r1,q2,r2,p2,k);
    return internal::intersection_test_vertex(p1,q1,r1,q2,r2,p2,k);
  }
  return internal::intersection_test_vertex(p1,q1,r1,r2,p2,q2,k);
  
}

} // namespace internal




template <class K>
inline bool do_intersect(const Triangle_2<K> &t1, 
			 const Triangle_2<K> &t2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(t1,t2);
}

} //namespace CGAL

#endif //CGAL_TRIANGLE_2_TRIANGLE_2_DO_INTERSECT_H

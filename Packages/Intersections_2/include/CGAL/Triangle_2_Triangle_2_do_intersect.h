// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Triangle_2_Triangle_2_do_intersect.h
// package       : 
// maintainer    : 
// source        : 
// author(s)     : Philippe Guigue
//
// coordinator   : 
//
// ======================================================================

#ifndef CGAL_TRIANGLE_2_TRIANGLE_2_DO_INTERSECT_H
#define CGAL_TRIANGLE_2_TRIANGLE_2_DO_INTERSECT_H


CGAL_BEGIN_NAMESPACE


template <class K>
bool do_intersect(const Triangle_2<K> &t1, 
		  const Triangle_2<K> &t2,
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
      return CGALi::intersection_test_edge(p1,q1,r1,p2,q2,r2,k);
    } 
    if ( orientation(*r2,*p2,*p1) != NEGATIVE ) 
      return CGALi::intersection_test_edge(p1,q1,r1,r2,p2,q2,k);
    return CGALi::intersection_test_vertex(p1,q1,r1,p2,q2,r2,k);
    
  }
  
  if ( orientation(*q2,*r2,*p1) != NEGATIVE ) {
    if ( orientation(*r2,*p2,*p1) != NEGATIVE ) 
      return CGALi::intersection_test_edge(p1,q1,r1,q2,r2,p2,k);
    return CGALi::intersection_test_vertex(p1,q1,r1,q2,r2,p2,k);
  }
  return CGALi::intersection_test_vertex(p1,q1,r1,r2,p2,q2,k);
  
}


namespace CGALi {

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

}

template <class K>
inline bool do_intersect(const Triangle_2<K> &t1, 
			 const Triangle_2<K> &t2){

  return do_intersect(t1,t2,K());
}

CGAL_END_NAMESPACE


#endif //CGAL_TRIANGLE_2_TRIANGLE_2_DO_INTERSECT_H





// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Triangle_3_Point_3_do_intersect.h
// package       : Intersections_3
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Philippe Guigue
//
// maintainer    : 
//
//
// ======================================================================

#ifndef CGAL_TRIANGLE_3_POINT_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_POINT_3_DO_INTERSECT_H


CGAL_BEGIN_NAMESPACE


template <class K>
bool do_intersect(const Triangle_3<K> &t, 
		  const Point_3<K>    &p,
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
inline bool do_intersect(const Point_3<K> &p, 
			 const Triangle_3<K> &t)
{ 
  return do_intersect(t,p,K());
} 

template <class K>
inline bool do_intersect(const Triangle_3<K> &t, 
			 const Point_3<K> &p)
{ 
  return do_intersect(t,p,K());
} 


template <class K>
inline bool do_intersect(const Point_3<K> &p, 
			 const Triangle_3<K> &t,
			 const K & k)
{ 
  return do_intersect(t,p,k);
} 

CGAL_END_NAMESPACE


#endif // CGAL_TRIANGLE_3_POINT_3_DO_INTERSECT_H

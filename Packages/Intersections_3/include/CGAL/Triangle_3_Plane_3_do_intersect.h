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
// file          : include/CGAL/Triangle_3_Plane_3_do_intersect.h
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

#ifndef CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H


CGAL_BEGIN_NAMESPACE


template <class K>
bool do_intersect(const Triangle_3<K> &t, 
		  const Plane_3<K>   &h,
		  const K & k)
{
  
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t)) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(h)) ;

   
  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();
  
  typename K::Construct_point_on_3 point_on = 
    k.construct_point_on_3_object();
  
  typename K::Oriented_side_3 oriented_side =
    k.oriented_side_3_object();

  

  switch ( oriented_side(h,vertex_on(t,0)) ) {
  case ON_POSITIVE_SIDE:
    return oriented_side(h,vertex_on(t,1)) != ON_POSITIVE_SIDE
      || oriented_side(h,vertex_on(t,2)) != ON_POSITIVE_SIDE;
  case ON_NEGATIVE_SIDE:
    return oriented_side(h,vertex_on(t,1)) != ON_NEGATIVE_SIDE
      || oriented_side(h,vertex_on(t,2)) != ON_NEGATIVE_SIDE ;
  case ON_ORIENTED_BOUNDARY:
    return true;
  default:// should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}

template <class K>
inline bool do_intersect(const Triangle_3<K> &t, 
			 const Plane_3<K>    &h)
{
  return do_intersect(t,h,K());
}
  
template <class K>
inline bool do_intersect(const Plane_3<K>    &h, 
			 const Triangle_3<K> &t)
{
  return do_intersect(t,h,K());
}

template <class K>
inline bool do_intersect(const Plane_3<K>    &h, 
			 const Triangle_3<K> &t,
			 const K & k)
{
  return do_intersect(t,h,k);
}



CGAL_END_NAMESPACE


#endif //CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H





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
// file          : include/CGAL/Triangle_3_Line_3_do_intersect.h
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

#ifndef CGAL_TRIANGLE_3_LINE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_LINE_3_DO_INTERSECT_H


CGAL_BEGIN_NAMESPACE


template <class K>
bool do_intersect(const Triangle_3<K> &t, 
		  const Line_3<K>     &l,
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
      case COLLINEAR:
	switch ( pqbc ) {
	case POSITIVE: return  orientation(p,q,c,a) != NEGATIVE ;
	case NEGATIVE: return  orientation(p,q,c,a) != POSITIVE ;
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
  
  // Coplanar case
    
  const Orientation pqa = coplanar_orientation(p,q,a);
  
  return  coplanar_orientation(p,q,b) != pqa 
    ||  coplanar_orientation(p,q,c) != pqa ;
  
}


template <class K>
bool do_intersect(const Triangle_3<K> &t, 
		  const Line_3<K>     &l)
{
  return do_intersect(t,l,K());
}


template <class K>
inline bool do_intersect(const Line_3<K> &l, 
			 const Triangle_3<K> &t)
{
  return do_intersect(t,l,K());
}


template <class K>
inline bool do_intersect(const Line_3<K> &l, 
			 const Triangle_3<K> &t,
			 const K  & k )
{
  return do_intersect(t,l,k);
}

CGAL_END_NAMESPACE


#endif //CGAL_TRIANGLE_3_LINE_3_DO_INTERSECT_H





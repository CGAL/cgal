// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/convexity_check_3.h
// package       : $CGAL_Package: Convex_hull_3 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Convex Hulls and Extreme Points
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//               : Amol Prakash <prakash@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: functions for testing validity of convex hulls
// ============================================================================

#ifndef CGAL_CONVEXITY_CHECK_3_H
#define CGAL_CONVEXITY_CHECK_3_H

#include <CGAL/intersections.h>

namespace CGAL {

template <class Facet_handle, class Traits>
bool is_locally_convex(Facet_handle f_hdl, const Traits& traits)
{
  // This function checks if all the faces around facet *f_hdl form part of 
  // the convex hull

  typedef typename Facet_handle::value_type                Facet;
  typedef typename Facet::Halfedge_around_facet_circulator Halfedge_circ;
  typedef typename Traits::Point_3			   Point_3; 

  typename Traits::Has_on_positive_side_3 has_on_positive_side = 
            traits.has_on_positive_side_3_object();
  
  Halfedge_circ h_circ = (*f_hdl).facet_begin();
 
  do 
  {
     // Take the point on the other facet not shared by this facet
     Point_3 point= h_circ->opposite()->next()->vertex()->point();
     
     // Point must be on the plane or on the negative side 
     if (has_on_positive_side((*f_hdl).plane(), point))
          return false;
     h_circ++;
  }
  while ( h_circ != (*f_hdl).facet_begin());
  return true;
}


// Pre: equations of facet planes have been computed
template<class Polyhedron, class Traits>
bool is_strongly_convex_3(Polyhedron& P, const Traits& traits)
{ 
  typedef typename Polyhedron::Facet_iterator  Facet_iterator;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Traits::Point_3             Point_3; 
  typedef typename Traits::Ray_3               Ray_3; 
  typedef typename Traits::Triangle_3          Triangle_3; 
  
  Facet_iterator f_it;
  for ( f_it = P.facets_begin(); f_it != P.facets_end(); ++f_it)
     if (!is_locally_convex(f_it, traits)) 
	return false;
 
  // Check 2: see if a point interior to the hull is actually on the same
  // side of each facet of P
 
  typename Traits::Coplanar_3 coplanar = traits.coplanar_3_object();

  Vertex_iterator v_it = P.vertices_begin();
  Point_3 p;
  Point_3 q;
  Point_3 r;
  Point_3 s;

  // First take 3 arbitrary points
  p = v_it->point();  v_it++;
  q = v_it->point();  v_it++;
  r = v_it->point();  v_it++;
  
  // Now take 4th point s.t. it's not coplaner with them
  do
  {
    s = v_it->point();  
    v_it++;
  } while (coplanar(p,q,r,s) && (v_it != P.vertices_end()));
  

  // if no such point, all are coplanar so it is not strongly convex
  if( v_it == P.vertices_end() )
    return false;

  
  // else construct a point inside the polyhedron
  typename Traits::Construct_centroid_3 construct_centroid =
           traits.construct_centroid_3_object();
  Point_3 inside_pt = construct_centroid(p,q,r,s);


  typename Traits::Oriented_side_3 oriented_side = 
            traits.oriented_side_3_object();

  f_it = P.facets_begin();

  Oriented_side side = oriented_side(f_it->plane(), inside_pt);

  // the point inside should not be on the facet plane
  if (side == ON_ORIENTED_BOUNDARY) return false;

  // now make sure this point that is inside the polyhedron is on the same
  // side of each facet
  for (f_it++; f_it != P.facets_end(); f_it++)
  {
    if ( oriented_side(f_it->plane(), inside_pt) != side )
       return false;
  }

  // Check 3 :  see if a ray from the interior point to a point in the
  // middle of one of the facets intersects any other facets 
  typename Traits::Construct_ray_3 construct_ray = 
            traits.construct_ray_3_object();
  typename Traits::Construct_triangle_3 construct_triangle = 
            traits.construct_triangle_3_object();
  typename Traits::Intersect_3 intersect =
            traits.intersect_3_object();
  typename Traits::Has_on_3 has_on =
            traits.has_on_3_object();
  f_it = P.facets_begin();
  Point_3 facet_pt = 
     construct_centroid(f_it->halfedge()->opposite()->vertex()->point(),
                        f_it->halfedge()->vertex()->point(),
                        f_it->halfedge()->next()->vertex()->point());
  Ray_3  ray = construct_ray(inside_pt, facet_pt);
  Point_3 i_point;
  Point_3 i_ray;
  for ( ++f_it ; f_it != P.facets_end(); f_it++)
  {
     Object obj = intersection(f_it->plane(), ray);
     if (assign(i_point, obj))  // ray  may intersect facet plane but not facet
     {
        Triangle_3 facet_tri = 
          construct_triangle(f_it->halfedge()->opposite()->vertex()->point(),
                             f_it->halfedge()->vertex()->point(),
                             f_it->halfedge()->next()->vertex()->point());
     
        if (has_on(facet_tri, i_point))
          return false;
     }
     else if (assign(i_ray, obj)) // ray should not lie on a facet plane
        return false;
  }
  return true;
}

template<class Polyhedron, class R>
bool CGAL_is_strongly_convex_3(Polyhedron& P, Point_3<R>*)
{ 
  return is_strongly_convex_3(P, R());
}

template<class Polyhedron>
bool is_strongly_convex_3(Polyhedron& P)
{ 
  typedef typename Polyhedron::Point_3 Point_3;

  return CGAL_is_strongly_convex_3(P, reinterpret_cast<Point_3*>(0));
}
  
template <class ForwardIterator, class Polyhedron, class Traits>
bool all_points_inside( ForwardIterator first, 
			ForwardIterator last,
			Polyhedron& P,
			const Traits&  traits)
{
   typedef typename Polyhedron::Facet_iterator Facet_iterator;
   typename Traits::Has_on_positive_side_3   has_on_positive_side = 
             traits.has_on_positive_side_3_object();

   for (ForwardIterator p_it = first; p_it != last; p_it++)
   {
      Facet_iterator f_it;
      for (f_it = P.facets_begin(); f_it != P.facets_end(); f_it++)
      {
	 if (has_on_positive_side(f_it->plane(),*p_it))
             return false;
      }
   }
   return true;
}

} // namespace CGAL

#endif // CGAL_CONVEXITY_CHECK_3_H

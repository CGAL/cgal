// Copyright (c) 2001  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//               : Amol Prakash <prakash@mpi-sb.mpg.de>

#ifndef CGAL_CONVEXITY_CHECK_3_H
#define CGAL_CONVEXITY_CHECK_3_H

#include <CGAL/license/Convex_hull_3.h>


#include <CGAL/intersections.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/property_maps.h>

namespace CGAL {


  template <class Polyhedron, class Vpmap, class Plane, class Facet_handle>
  void get_plane2(const Polyhedron& P, Vpmap vpmap, Plane& plane, Facet_handle f)
{
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor h = halfedge(f,P);
  plane = Plane(get(vpmap, target(opposite(h,P),P)),
                get(vpmap, target(h,P)),
                get(vpmap, target(next(h,P),P)));
}


  template <class Polyhedron, class Vpmap, class Facet_handle, class Traits>
bool is_locally_convex(Polyhedron P, Vpmap vpmap, Facet_handle f_hdl, const Traits& traits)
{
  // This function checks if all the faces around facet *f_hdl form part of
  // the convex hull

  typedef Halfedge_around_face_circulator<Polyhedron>      Halfedge_circ;
  typedef typename Traits::Point_3                           Point_3;
  typedef typename Traits::Plane_3                           Plane_3;

  typename Traits::Has_on_positive_side_3 has_on_positive_side =
            traits.has_on_positive_side_3_object();

  Halfedge_circ h_circ( halfedge(f_hdl,P),P), done(h_circ);

  do
  {
     // Take the point on the other facet not shared by this facet
    Point_3 point = get(vpmap, target(next(opposite(*h_circ,P),P),P));
     Plane_3 plane;
     get_plane2(P, vpmap, plane, f_hdl);
     // Point must be on the plane or on the negative side
     if (has_on_positive_side(plane, point)) {
       return false;
     }
     h_circ++;
  }
  while ( h_circ != done);
  return true;
}


// Pre: equations of facet planes have been computed
template<class Polyhedron, class Traits>
bool is_strongly_convex_3(const Polyhedron& P, const Traits& traits)
{
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::face_iterator face_iterator;

  typedef typename Traits::Point_3             Point_3;
  typedef typename Traits::Ray_3               Ray_3;
  typedef typename Traits::Triangle_3          Triangle_3;
  typedef typename Traits::Plane_3             Plane_3;

  typename boost::property_map<Polyhedron, vertex_point_t>::const_type vpmap  = get(CGAL::vertex_point, P);

  vertex_iterator v_it, v_it_e;
  boost::tie(v_it, v_it_e) = vertices(P);

  if (v_it == v_it_e) return false;

  for(face_descriptor fd : faces(P))
    if (!is_locally_convex(P, vpmap, fd, traits))
        return false;

  // Check 2: see if a point interior to the hull is actually on the same
  // side of each facet of P

  typename Traits::Coplanar_3 coplanar = traits.coplanar_3_object();

  face_iterator f_it, f_it_e;
  boost::tie(f_it, f_it_e) = faces(P);
  Point_3 p;
  Point_3 q;
  Point_3 r;
  Point_3 s;

  // First take 3 arbitrary points
  p = get(vpmap, *v_it);  v_it++;
  q = get(vpmap,*v_it);  v_it++;
  r = get(vpmap,*v_it);  v_it++;

  // three vertices that form a single (triangular) facet
  if (v_it == v_it_e){
    return f_it != f_it_e;
  }

  // Now take 4th point s.t. it's not coplaner with them
  while (v_it != v_it_e && coplanar(p, q, r, get(vpmap,*v_it)))
    v_it++;

  // if no such point, all are coplanar so it is not strongly convex
  if( v_it == v_it_e ){
    return false;
  }

  s = get(vpmap,*v_it);

  // else construct a point inside the polyhedron
  typename Traits::Construct_centroid_3 construct_centroid =
           traits.construct_centroid_3_object();
  Point_3 inside_pt = construct_centroid(p,q,r,s);

  typename Traits::Oriented_side_3 oriented_side =
            traits.oriented_side_3_object();


  Plane_3 plane;
  get_plane2(P, vpmap, plane, *f_it);
  Oriented_side side = oriented_side(plane, inside_pt);

  // the point inside should not be on the facet plane
  if (side == ON_ORIENTED_BOUNDARY){
    return false;
  }

  // now make sure this point that is inside the polyhedron is on the same
  // side of each facet
  for (f_it++; f_it != f_it_e; f_it++)
  {
    Plane_3 plane;
    get_plane2(P, vpmap, plane, *f_it);
    if ( oriented_side(plane, inside_pt) != side ){
      return false;
    }
  }


  // Check 3 :  see if a ray from the interior point to a point in the
  // middle of one of the facets intersects any other facets
  typename Traits::Construct_ray_3 construct_ray =
            traits.construct_ray_3_object();
  typename Traits::Construct_triangle_3 construct_triangle =
            traits.construct_triangle_3_object();
  typename Traits::Do_intersect_3 do_intersect =
            traits.do_intersect_3_object();

  f_it = faces(P).first;

  Point_3 facet_pt =
    construct_centroid(get(vpmap,target(opposite(halfedge(*f_it,P),P),P)),
                       get(vpmap,target(halfedge(*f_it,P),P)),
                       get(vpmap,target(next(halfedge(*f_it,P),P),P)));
  Ray_3  ray = construct_ray(inside_pt, facet_pt);

  for ( ++f_it ; f_it != f_it_e; f_it++)
  {
    Triangle_3 facet_tri =
      construct_triangle(get(vpmap,target(opposite(halfedge(*f_it,P),P),P)),
                         get(vpmap,target(halfedge(*f_it,P),P)),
                         get(vpmap,target(next(halfedge(*f_it,P),P),P)));

    if (do_intersect(facet_tri, ray)){
      return false;
    }
  }

  return true;
}

template<class Polyhedron, class R>
bool CGAL_is_strongly_convex_3(const Polyhedron& P, Point_3<R>*)
{
  return is_strongly_convex_3(P, R());
}

template<class Polyhedron>
bool is_strongly_convex_3(const Polyhedron& P)
{
  typedef typename boost::property_map<Polyhedron, vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point_3;

  return CGAL_is_strongly_convex_3(P, reinterpret_cast<Point_3*>(0));
}

template <class ForwardIterator, class Polyhedron, class Traits>
bool all_points_inside( ForwardIterator first,
                        ForwardIterator last,
                        const Polyhedron& P,
                        const Traits&  traits)
{
  typedef typename Traits::Plane_3                           Plane_3;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
   typename Traits::Has_on_positive_side_3   has_on_positive_side =
             traits.has_on_positive_side_3_object();

   for (ForwardIterator p_it = first; p_it != last; p_it++)
   {
      for(face_descriptor fd : faces(P))
      {
        Plane_3 plane;
        get_plane2(P,plane, fd);
        if (has_on_positive_side(plane,*p_it)){
             return false;
        }
      }
   }
   return true;
}

} // namespace CGAL

#endif // CGAL_CONVEXITY_CHECK_3_H

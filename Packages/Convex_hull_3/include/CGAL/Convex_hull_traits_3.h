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
// file          : include/CGAL/Convex_hull_traits_3.h
// package       : $CGAL_Package: Convex_hull_3 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Convex Hulls and Extreme Points
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: 3D convex hull traits class
// ============================================================================

#ifndef CGAL_CONVEX_HULL_TRAITS_3_H
#define CGAL_CONVEX_HULL_TRAITS_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/distance_predicates_3.h>
#include <CGAL/ch_predicate_classes_3.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_projective_xy_traits_2.h>
#include <CGAL/Convex_hull_projective_xz_traits_2.h>
#include <CGAL/Convex_hull_projective_yz_traits_2.h>


namespace CGAL {

/*
template <class R_>
class Convex_hull_traits_3: public R_
{
 public:  
  typedef   R_                                              R;
  typedef   typename R::Point_3                             Point_3;
  typedef   typename R::Plane_3                             Plane_3;
  typedef   Less_signed_distance_to_plane_3<Plane_3, Point_3>
                                                Less_signed_distance_to_plane_3;
  typedef   Construct_centroid_3<Point_3> Construct_centroid_3;

  Construct_centroid_3  
  construct_centroid_3_object() const
  { return Construct_centroid_3(); }

  Less_signed_distance_to_plane_3  
  less_signed_distance_to_plane_3_object(Plane_3& P) const
  { return Less_signed_distance_to_plane_3(P); }

};

*/

template <class R_>
class Convex_hull_traits_3 
{
 public:  
  typedef R_                                     R;
  typedef typename R::Point_3                    Point_3;
  typedef typename R::Segment_3                  Segment_3;
  typedef typename R::Triangle_3                 Triangle_3;
  typedef typename R::Plane_3                    Plane_3;
  typedef typename R::Vector_3                   Vector_3;

  typedef Polyhedron_default_traits_3<R>         Polyhedron_traits;
  typedef Halfedge_data_structure_polyhedron_default_3<R>         
                                                  HDS;
  typedef Polyhedron_3<Polyhedron_traits, HDS>   Polyhedron_3;

  typedef typename R::Construct_segment_3        Construct_segment_3;
  typedef typename R::Construct_ray_3            Construct_ray_3;
  typedef typename R::Construct_plane_3          Construct_plane_3;
  typedef typename R::Construct_vector_3         Construct_vector_3;
  typedef typename R::Construct_triangle_3       Construct_triangle_3;
  typedef Construct_centroid_3<Point_3>          Construct_centroid_3;
  typedef Construct_orthogonal_vector_3<Plane_3, Vector_3> 
                                                 Construct_orthogonal_vector_3;

  typedef typename R::Collinear_3                Collinear_3;
  typedef typename R::Coplanar_3                 Coplanar_3;
  typedef typename R::Less_distance_to_point_3   Less_distance_to_point_3;
  typedef typename R::Has_on_positive_side_3     Has_on_positive_side_3;

  typedef Less_signed_dist_to_plane_3<Plane_3, Point_3>
                                               Less_signed_distance_to_plane_3;

  // required for degenerate case of all points coplanar
  typedef Convex_hull_projective_xy_traits_2<Point_3>  Traits_xy;
  typedef Convex_hull_projective_xz_traits_2<Point_3>  Traits_xz;
  typedef Convex_hull_projective_yz_traits_2<Point_3>  Traits_yz;

  // for postcondition checking 
  typedef typename R::Ray_3                      Ray_3; 

  typedef typename R::Has_on_3                   Has_on_3;
  typedef typename R::Oriented_side_3            Oriented_side_3;
  typedef typename R::Intersect_3                Intersect_3;

  Construct_segment_3
  construct_segment_3_object() const
  { return Construct_segment_3(); }

  Construct_ray_3
  construct_ray_3_object() const
  { return Construct_ray_3(); }

  Construct_plane_3
  construct_plane_3_object() const
  { return Construct_plane_3(); }

  Construct_triangle_3
  construct_triangle_3_object() const
  { return Construct_triangle_3(); }

  Construct_vector_3
  construct_vector_3_object() const
  { return Construct_vector_3(); }

  Construct_centroid_3  
  construct_centroid_3_object() const
  { return Construct_centroid_3(); }

  Construct_orthogonal_vector_3
  construct_orthogonal_vector_3_object() const
  { return Construct_orthogonal_vector_3(); }

  Collinear_3
  collinear_3_object() const
  { return Collinear_3(); }

  Coplanar_3
  coplanar_3_object() const
  { return Coplanar_3(); }
 
  Has_on_3
  has_on_3_object() const
  { return Has_on_3(); }

  Less_distance_to_point_3
  less_distance_to_point_3_object(const Point_3& p) const
  { return Less_distance_to_point_3(p); }

  Has_on_positive_side_3
  has_on_positive_side_3_object() const
  { return Has_on_positive_side_3(); }

  Oriented_side_3
  oriented_side_3_object() const
  { return Oriented_side_3(); }

  Intersect_3
  intersect_3_object() const
  { return Intersect_3(); }

  Less_signed_distance_to_plane_3  
  less_signed_distance_to_plane_3_object(Plane_3 p) const
  { return Less_signed_distance_to_plane_3(p); }
};

} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_3_H



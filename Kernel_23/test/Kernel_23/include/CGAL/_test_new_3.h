// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel
//                 Stefan Schirra
 

#ifndef CGAL_TEST_NEW_3_H
#define CGAL_TEST_NEW_3_H

#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/_test_compare_dihedral_angle_3.h>

#include <CGAL/Testsuite/use.h>

using CGAL::internal::use;

// Accessory function testing functions that require sqrt().
// Doesn't instantiate anything if RT doesn't support sqrt().
template <class R>
bool
_test_new_3_sqrt(const R&, CGAL::Tag_false)
{
    //bool UNTESTED_STUFF_BECAUSE_SQRT_IS_NOT_SUPPORTED;
    std::cout << std::endl
              << "NOTE : FT doesn't support sqrt(),"
                 " hence some functions are not tested."
	      << std::endl;
    return true;
}

template <class R>
bool
_test_new_3_sqrt(const R& rep, CGAL::Tag_true)
{
  typedef typename R::FT          FT;
  typedef typename R::Point_3     Point_3;
  typedef typename R::Vector_3    Vector_3;
  typedef typename R::Triangle_3  Triangle_3;

  Point_3 p3 (1,1,1);
  Point_3 p4 (1,1,2,2);
  Point_3 p5 (1,2,3,4);

  typename R::Construct_triangle_3 construct_triangle
        = rep.construct_triangle_3_object();
  Triangle_3 t2 = construct_triangle(p3,p4,p5);

  typename R::Compute_area_3 compute_area
        = rep.compute_area_3_object();
  FT tmp11a = compute_area(t2);
     tmp11a = compute_area(p3, p4, p5);


  typename R::Construct_unit_normal_3 construct_unit_normal
    = rep.construct_unit_normal_3_object();

  Vector_3 tmp11b = construct_unit_normal(Point_3(CGAL::ORIGIN), p4, p5);

  use(tmp11a);
  use(tmp11b);

  return true;
}


template <class R>
bool
test_new_3(const R& rep)
{
  std::cout << "Testing 3 dimensional functionality" << std::endl;

  using namespace CGAL;

  typedef typename R::RT                          RT;
  typedef typename R::FT                          FT;

  typedef typename R::Point_3                     Point_3;
  typedef typename R::Vector_3                    Vector_3;
  typedef typename R::Direction_3                 Direction_3;
  typedef typename R::Segment_3                   Segment_3;
  typedef typename R::Line_3                      Line_3;
  typedef typename R::Ray_3                       Ray_3;
  typedef typename R::Plane_3                     Plane_3;
  typedef typename R::Sphere_3                    Sphere_3;
  typedef typename R::Triangle_3                  Triangle_3;
  typedef typename R::Tetrahedron_3               Tetrahedron_3;
  typedef typename R::Iso_cuboid_3                Iso_cuboid_3;
  typedef typename R::Object_3                    Object_3;
  typedef typename R::Point_2                     Point_2;

  typename R::Construct_point_3 construct_point
        = rep.construct_point_3_object();
  Point_3 p1;
  Point_3 p2 = construct_point(ORIGIN);
  Point_3 p3 = construct_point(1,1,1);
  Point_3 p3bis = construct_point(RT(1),RT(1),RT(1));
  Point_3 p3ter = construct_point(FT(1),FT(1),FT(1));
  use(p3bis); use(p3ter);
  Point_3 p4 = construct_point(1,1,2,2);
  Point_3 p5 = construct_point(1,2,3,4);
  Point_3 p6 = construct_point(4,2,1,2);

  typename R::Construct_vector_3 construct_vector
        = rep.construct_vector_3_object();
  Vector_3 v1;
  Vector_3 v2 = construct_vector(NULL_VECTOR);
  Vector_3 v3 = construct_vector(1,1,1);
  Vector_3 v3bis = construct_vector(RT(1),RT(1),RT(1));
  Vector_3 v3ter = construct_vector(FT(1),FT(1),FT(1));
  use(v3bis); use(v3ter);
  Vector_3 v4 = construct_vector(1,1,2,2);
  Vector_3 v5 = construct_vector(p5, p6);

  typename R::Construct_direction_3 construct_direction
        = rep.construct_direction_3_object();
  Direction_3 d1;
  Direction_3 d2 = construct_direction(v3);
  Direction_3 d3 = construct_direction(1,1,5);
  Direction_3 d4 = construct_direction(1,5,5);
  // remaining constructions tested below, after the 
  // corresponding types have been introduced

  typename R::Construct_segment_3 construct_segment
        = rep.construct_segment_3_object();
  Segment_3 s1;
  Segment_3 s2 = construct_segment(p2,p3);

  typename R::Construct_ray_3 construct_ray =
        rep.construct_ray_3_object();
  Ray_3 r1;
  Ray_3 r2 = construct_ray(p2,p4);
  Ray_3 r3 = construct_ray(p2,d3);
  Ray_3 r4 = construct_ray(p2,v3);

  typename R::Construct_line_3 construct_line
        = rep.construct_line_3_object();
  Line_3 l1;
  Line_3 l2 = construct_line(p5,p6);
  Line_3 l3 = construct_line(p2,p3);
  Line_3 l4 = construct_line(p2,d4);
  Line_3 l5 = construct_line(s2);
  Line_3 l6 = construct_line(r2);
  Line_3 l7 = construct_line(p2, construct_direction(v4));
  Line_3 l8 = construct_line(p2,v4);

  // remaining construct_direction tests
  Direction_3 d5 = construct_direction(l3);
  Direction_3 d6 = construct_direction(r2);
  Direction_3 d7 = construct_direction(s2);

  // remaining construct_ray tests
  Ray_3 r5 = construct_ray(p2, l3);

  // remaining construct_vector tests
  Vector_3 v7 = construct_vector(s2);
  Vector_3 v8 = construct_vector(r2);
  Vector_3 v9 = construct_vector(l2);
  
  typename R::Construct_plane_3 construct_plane
        = rep.construct_plane_3_object();
  Plane_3 h1;
  Plane_3 h2 = construct_plane(1,1,1,1);
  Plane_3 h3 = construct_plane(p2,p3,p4);
  Plane_3 h4 = construct_plane(p2,d4);
  Plane_3 h5 = construct_plane(l2,p4);
  Plane_3 h6 = construct_plane(r2,p4);
  Plane_3 h7 = construct_plane(s2,p4);
  Plane_3 h8 = construct_plane(p2,v3);

  typename R::Construct_sphere_3 construct_sphere
        = rep.construct_sphere_3_object();
  Sphere_3 sp1 = construct_sphere(p2,1);
  Sphere_3 sp2 = construct_sphere(p2,1,COUNTERCLOCKWISE);
  Sphere_3 sp3 = construct_sphere(p2,p3,p4,p5);
  Sphere_3 sp4 = construct_sphere(p2,p3,p4);
  Sphere_3 sp5 = construct_sphere(p2,p3,p4,CLOCKWISE);
  Sphere_3 sp6 = construct_sphere(p2,p3);
  Sphere_3 sp7 = construct_sphere(p2,p3,CLOCKWISE);
  Sphere_3 sp8 = construct_sphere(p3);
  Sphere_3 sp9 = construct_sphere(p3,CLOCKWISE);


  typename R::Construct_triangle_3 construct_triangle
        = rep.construct_triangle_3_object();
  Triangle_3 t1;
  Triangle_3 t2 = construct_triangle(p2,p3,p4);

  typename R::Construct_tetrahedron_3 construct_tetrahedron
        = rep.construct_tetrahedron_3_object();
  Tetrahedron_3 th1;
  Tetrahedron_3 th2 = construct_tetrahedron(p2,p3,p4,p5);

  typename R::Construct_iso_cuboid_3 construct_iso_cuboid
        = rep.construct_iso_cuboid_3_object();
  Iso_cuboid_3 iso1 = construct_iso_cuboid(p3,p6);
               iso1 = construct_iso_cuboid(p2,p3,0);
               iso1 = construct_iso_cuboid(p3,p3,p6,p6,p4,p4);

  typename R::Construct_object_3 construct_object
        = rep.construct_object_3_object();
  Object_3 obj = construct_object(iso1);
           obj = construct_object(th2);
           obj = construct_object(t2);

  typename R::Construct_point_on_3 construct_point_on
        = rep.construct_point_on_3_object();
  Point_3 tmp1 = construct_point_on(l2, 0);

  typename R::Construct_projected_point_3 construct_projected_point
        = rep.construct_projected_point_3_object();
          tmp1 = construct_projected_point(l2, p4);
          tmp1 = construct_projected_point(h7, p3);

  typename R::Construct_lifted_point_3 construct_lifted_point
        = rep.construct_lifted_point_3_object();
          tmp1 = construct_lifted_point(h7, Point_2(10,4));

  typename R::Construct_scaled_vector_3 construct_scaled_vector
        = rep.construct_scaled_vector_3_object();
  Vector_3 v6 = construct_scaled_vector(v5, RT(5));
           v6 = construct_scaled_vector(v5, FT(5));

  typename R::Construct_translated_point_3 construct_translated_point
        = rep.construct_translated_point_3_object();
          p1 = construct_translated_point(tmp1, v6);
          p2 = construct_translated_point(p3, -v6);

  typename R::Construct_vertex_3 construct_vertex_3
        = rep.construct_vertex_3_object();
  Point_3 tmp2f = construct_vertex_3(s2, 0);
          tmp2f = construct_vertex_3(iso1, 0);
          tmp2f = construct_vertex_3(t2, 0);
          tmp2f = construct_vertex_3(th2, 0);
  
  typename R::Construct_min_vertex_3 construct_min_vertex_3
        = rep.construct_min_vertex_3_object();
          tmp2f = construct_min_vertex_3(iso1);
          tmp2f = construct_min_vertex_3(s2);

  typename R::Construct_max_vertex_3 construct_max_vertex_3
        = rep.construct_max_vertex_3_object();
          tmp2f = construct_max_vertex_3(iso1);
          tmp2f = construct_max_vertex_3(s2);


  typename R::Construct_normal_3 construct_normal
    = rep.construct_normal_3_object();

  Vector_3 tmp2g = construct_normal(Point_3(CGAL::ORIGIN), Point_3(1,0,0), Point_3(0,1,0));
  use(tmp2g);
  typename R::Construct_bbox_3 construct_bbox_3
    = rep.construct_bbox_3_object();

  Bbox_3 bb1 = construct_bbox_3(p1); // Point_3
  Bbox_3 bb2 = construct_bbox_3(s1); // Segment_3
  Bbox_3 bb3 = construct_bbox_3(t1); // Triangle_3
  Bbox_3 bb4 = construct_bbox_3(th1); // Tetrahedron_3
  Bbox_3 bb5 = construct_bbox_3(sp1); // Sphere_3
  Bbox_3 bb6 = construct_bbox_3(iso1); // Iso_cuboid_3

  typename R::Construct_cartesian_const_iterator_3 
    construct_cartesian_const_iterator_3
    = rep.construct_cartesian_const_iterator_3_object();

  typename R::Cartesian_const_iterator_3 cccit;

  cccit = construct_cartesian_const_iterator_3(p1);
  cccit = construct_cartesian_const_iterator_3(p1,0);
  cccit = construct_cartesian_const_iterator_3(v5);
  cccit = construct_cartesian_const_iterator_3(v5,0);

  typename R::Construct_perpendicular_plane_3 construct_perpendicular_plane
        = rep.construct_perpendicular_plane_3_object();
  Plane_3 tmp3 = construct_perpendicular_plane(l2,p2);

  typename R::Construct_equidistant_line_3 construct_equidistant_line
        = rep.construct_equidistant_line_3_object();
  Line_3 tmp3a = construct_equidistant_line(p1, p2, p3);

  typename R::Construct_perpendicular_line_3 construct_perpendicular_line
        = rep.construct_perpendicular_line_3_object();
         tmp3a = construct_perpendicular_line(h2,p5);

  typename R::Construct_orthogonal_vector_3 construct_orthogonal_vector
        = rep.construct_orthogonal_vector_3_object();
  Vector_3 tmp3b = construct_orthogonal_vector(h7);
  tmp3b =  construct_orthogonal_vector(p2, p3, p4);

  typename R::Construct_base_vector_3 construct_base_vector
        = rep.construct_base_vector_3_object();
           tmp3b = construct_base_vector(h7, 1);
           tmp3b = construct_base_vector(h7, 2);

  typename R::Construct_midpoint_3 construct_midpoint
        = rep.construct_midpoint_3_object();
  Point_3 tmp4 = construct_midpoint(p2,p3);

  typename R::Construct_center_3 construct_center
        = rep.construct_center_3_object();
          tmp4 = construct_center(sp2);

  typename R::Construct_circumcenter_3 construct_circumcenter
        = rep.construct_circumcenter_3_object();
          tmp4 = construct_circumcenter(p2,p3,p4,p5);
          tmp4 = construct_circumcenter(p2,p3,p4);
          tmp4 = construct_circumcenter(p2,p3);
          tmp4 = construct_circumcenter(th2);
          tmp4 = construct_circumcenter(t2);

  typename R::Construct_centroid_3 construct_centroid
        = rep.construct_centroid_3_object();
          tmp4 = construct_centroid(p2,p3,p4);
          tmp4 = construct_centroid(p2,p3,p4,p5);
          tmp4 = construct_centroid(t2);
          tmp4 = construct_centroid(th2);

  typename R::Construct_barycenter_3 construct_barycenter
        = rep.construct_barycenter_3_object();
          tmp4 = construct_barycenter(p2, FT(1), p3);
          tmp4 = construct_barycenter(p2, FT(1), p3, FT(2));
          tmp4 = construct_barycenter(p2, FT(1), p3, FT(2), p4);
          tmp4 = construct_barycenter(p2, FT(1), p3, FT(2), p4, FT(3));
          tmp4 = construct_barycenter(p2, FT(1), p3, FT(2), p4, FT(3), p5);
          tmp4 = construct_barycenter(p2, FT(1), p3, FT(2), p4, FT(3), p5, FT(4));

  typename R::Construct_cross_product_vector_3 construct_cross_product
        = rep.construct_cross_product_vector_3_object();
  Vector_3 tmp9 = construct_cross_product(v3,v4);


  typename R::Construct_opposite_direction_3 construct_opposite_direction
        = rep.construct_opposite_direction_3_object();
  Direction_3 tmp14a = construct_opposite_direction(d3);


  typename R::Construct_opposite_segment_3 construct_opposite_segment
        = rep.construct_opposite_segment_3_object();
  Segment_3 tmp5 = construct_opposite_segment(s2);


  typename R::Construct_opposite_ray_3 construct_opposite_ray
        = rep.construct_opposite_ray_3_object();
  Ray_3 tmp6 = construct_opposite_ray(r2);


  typename R::Construct_opposite_line_3 construct_opposite_line
        = rep.construct_opposite_line_3_object();
  Line_3 tmp7 = construct_opposite_line(l2);

  typename R::Construct_opposite_plane_3 construct_opposite_plane
        = rep.construct_opposite_plane_3_object();
  Plane_3 tmp71 = construct_opposite_plane(h2);

  typename R::Construct_opposite_sphere_3 construct_opposite_sphere
        = rep.construct_opposite_sphere_3_object();
  Sphere_3 sp1a = construct_opposite_sphere(sp1);

  typename R::Construct_opposite_vector_3 construct_opposite_vector
        = rep.construct_opposite_vector_3_object();
  Vector_3 tmp72 = construct_opposite_vector(v2);

  typename R::Construct_supporting_plane_3 construct_supporting_plane
        = rep.construct_supporting_plane_3_object();
  Plane_3 tmp8 = construct_supporting_plane(t2);


  typename R::Intersect_3 intersect
        = rep.intersect_3_object();
  Object_3 tmp10a = intersect(l2,h2);
  Object_3 tmp10b = intersect(r2,h2);

  bool tmp12a;
  bool tmp12b;

  typename R::Do_intersect_3 do_intersect
        = rep.do_intersect_3_object();
     tmp12a = do_intersect(l2,h2);
     tmp12b = do_intersect(t2,th2);
     tmp12b = do_intersect(th2,t2);
     tmp12b = do_intersect(r2,h2);

  typename R::Assign_3  assign
        = rep.assign_3_object();
       tmp12a = assign(p1,tmp10a);
       tmp12b = assign(p1,tmp10b);
  (void) tmp12a;
  (void) tmp12b;

  typename R::Compute_determinant_3 compute_determinant_3
        = rep.compute_determinant_3_object();
  FT tmp12aa = compute_determinant_3(v3, v4, v5);

  typename R::Compute_scalar_product_3 Compute_scalar_product
        = rep.compute_scalar_product_3_object();
  FT tmp12bb = Compute_scalar_product(v3, v4);
  (void) tmp12bb;

  typename R::Compute_squared_distance_3 Compute_squared_distance
        = rep.compute_squared_distance_3_object();
  FT tmp12c = Compute_squared_distance(p1, p2);
     tmp12c = Compute_squared_distance(p1, r2);
     tmp12c = Compute_squared_distance(p1, h2);
  (void) tmp12c;

  typename R::Compute_squared_length_3 compute_squared_length
        = rep.compute_squared_length_3_object();
  FT tmp11 = compute_squared_length(v3);
  tmp11 = compute_squared_length(s2);
  (void) tmp11;

  
  typename R::Compute_squared_radius_3 Compute_squared_radius
        = rep.compute_squared_radius_3_object();
  FT tmp11aa = Compute_squared_radius(sp1);
     tmp11aa = Compute_squared_radius(p3);
     tmp11aa = Compute_squared_radius(p3, p4);
     tmp11aa = Compute_squared_radius(p3, p4, p5);
     tmp11aa = Compute_squared_radius(p3, p4, p5, p6);
  (void) tmp11aa;

  typename R::Compute_squared_area_3 compute_squared_area
        = rep.compute_squared_area_3_object();
  FT tmp11a = compute_squared_area(t2);
     tmp11a = compute_squared_area(p3, p4, p5);
  (void) tmp11a;

  typename R::Compute_volume_3 compute_volume
        = rep.compute_volume_3_object();
  FT tmp11b = compute_volume(th2);
     tmp11b = compute_volume(iso1);
     tmp11b = compute_volume(p3, p4, p5, p6);
  (void) tmp11b;

  typename R::Equal_3 equal
        = rep.equal_3_object();
       tmp12a = equal(p2,p3);
       tmp12b = equal(l2,l3);
       bool tmp12d = equal(d2,d3);
       bool tmp12e = equal(s2,s2);
       bool tmp12f = equal(r2,r3);
       bool tmp12g = equal(h2,h3);
       bool tmp12h = equal(sp2,sp3);
       bool tmp12i = equal(t2,t2);
       bool tmp12j = equal(th2,th2);
       bool tmp12k = equal(iso1,iso1);
       bool tmp12l = equal(v2,v3);
  (void) tmp12a;
  (void) tmp12b;
  (void) tmp12d;
  (void) tmp12e;
  (void) tmp12f;
  (void) tmp12g;
  (void) tmp12h;
  (void) tmp12i;
  (void) tmp12j;
  (void) tmp12k;
  (void) tmp12l;

  typename R::Equal_x_3 equal_x
        = rep.equal_x_3_object();
  bool tmp13 = equal_x(p2,p3);
  (void) tmp13;


  typename R::Equal_y_3 equal_y
        = rep.equal_y_3_object();
  bool tmp14 = equal_y(p2,p3);
  (void) tmp14;


  typename R::Equal_z_3 equal_z
        = rep.equal_z_3_object();
  bool tmp15 = equal_z(p2,p3);
  (void) tmp15;


  typename R::Equal_xy_3 equal_xy
        = rep.equal_xy_3_object();
  bool tmp16 = equal_xy(p2,p3);
  (void) tmp16;

  typename R::Less_x_3 less_x
        = rep.less_x_3_object();
  bool tmp18 = less_x(p2,p3);
  (void) tmp18;

  typename R::Less_y_3 less_y
        = rep.less_y_3_object();
  bool tmp19 = less_y(p2,p3);
  (void) tmp19;

  typename R::Less_z_3 less_z
        = rep.less_z_3_object();
  bool tmp20 = less_z(p2,p3);
  (void) tmp20;

  typename R::Less_xy_3 less_xy
        = rep.less_xy_3_object();
  bool tmp21 = less_xy(p2,p3);
  (void) tmp21;

  typename R::Less_xyz_3 less_xyz
        = rep.less_xyz_3_object();
  bool tmp22 = less_xyz(p2,p3);
  (void) tmp22;

  typename R::Compare_x_3 compare_x
        = rep.compare_x_3_object();
  Comparison_result tmp23 = compare_x(p2,p3);
  (void) tmp23;


  typename R::Compare_y_3 compare_y
        = rep.compare_y_3_object();
  Comparison_result tmp24 = compare_y(p2,p3);
  (void) tmp24;


  typename R::Compare_z_3 compare_z
        = rep.compare_z_3_object();
  Comparison_result tmp25 = compare_z(p2,p3);
  (void) tmp25;


  typename R::Compare_xy_3 compare_xy
        = rep.compare_xy_3_object();
  Comparison_result tmp26 = compare_xy(p2,p3);
  (void) tmp26;


  typename R::Compare_xyz_3 compare_xyz
        = rep.compare_xyz_3_object();
  Comparison_result tmp27 = compare_xyz(p2,p3);
  (void) tmp27;

  typename R::Less_distance_to_point_3 less_distance_to_point
        = rep.less_distance_to_point_3_object();
  bool tmp28 = less_distance_to_point(p4,p2,p3);
  (void) tmp28;

  typename R::Less_signed_distance_to_plane_3 less_signed_distance_to_plane
        = rep.less_signed_distance_to_plane_3_object();
  bool tmp28a = less_signed_distance_to_plane(tmp8,p2,p3);
  (void) tmp28a;

  {
    bool tmp = _test_compare_dihedral_angle_3(rep);
    assert(tmp);
  }

  typename R::Compare_distance_3 compare_dist
        = rep.compare_distance_3_object();
  Comparison_result tmp34ab = compare_dist(p2,p3,p4);
  tmp34ab = compare_dist(p2,p3,p2,p3);
  tmp34ab = compare_dist(p1, p2, p3, p4);  
  (void) tmp34ab;

  typename R::Compare_squared_distance_3 compare_sq_dist
        = rep.compare_squared_distance_3_object();
  tmp34ab = compare_sq_dist(p2,p3,FT(1));
  tmp34ab = compare_sq_dist(p2,p3,p2,p3);
  tmp34ab = compare_sq_dist(p1, l2, FT(1));
  tmp34ab = compare_sq_dist(p2, p1, FT(1));
  tmp34ab = compare_sq_dist(l2, l3, FT(1));
  tmp34ab = compare_sq_dist(p1, s2, FT(1));
  tmp34ab = compare_sq_dist(p1, p2, p3, p4);

  tmp34ab = CGAL::compare_distance(p2,p3,p2,p3);
  tmp34ab = CGAL::compare_distance(p1, p2, p3, p4);  
  tmp34ab = CGAL::compare_distance(p1, p2, p3);  

  typename R::Compare_squared_radius_3 compare_sq_radius
        = rep.compare_squared_radius_3_object();

  {
    FT rad(0);
    Comparison_result tmp;
    tmp = compare_sq_radius(p1, p2, p3, p4, rad);
    tmp = compare_sq_radius(p1, p2, p3, rad);
    tmp = compare_sq_radius(p1, p2, rad);
    tmp = compare_sq_radius(p1, rad);
    (void)tmp;
  }

  typename R::Collinear_3 collinear
        = rep.collinear_3_object();
  bool tmp29 = collinear(p2,p3,p4);
  (void) tmp29;

  typename R::Coplanar_3 coplanar
        = rep.coplanar_3_object();
  bool tmp30 = coplanar(p2,p3,p4,p5);
  (void) tmp30;

  Point_3 p7(0,0,0);
  Point_3 p8(1,0,0);
  Point_3 p9(1,1,0);
  Point_3 p10(0,1,0);
  typename R::Coplanar_orientation_3 coplanar_orientation
        = rep.coplanar_orientation_3_object();
  Orientation tmp30a = coplanar_orientation(p7,p8,p9,p10);
              tmp30a = coplanar_orientation(p7,p8,p9);
  (void) tmp30a;

  typename R::Coplanar_side_of_bounded_circle_3
           coplanar_side_of_bounded_circle
        = rep.coplanar_side_of_bounded_circle_3_object();
  Bounded_side tmp30b = coplanar_side_of_bounded_circle(p7,p8,p9,p10);
  (void) tmp30b;

  typename R::Orientation_3 orientation
        = rep.orientation_3_object();
  Orientation tmp31 = orientation(p2,p3,p4,p5);
  tmp31 = orientation(v3,v4,v5);
  (void) tmp31;


  typename R::Is_degenerate_3 is_degenerate
        = rep.is_degenerate_3_object();
  bool tmp32 = is_degenerate(l2);
       tmp32 = is_degenerate(iso1);
       tmp32 = is_degenerate(h2);
       tmp32 = is_degenerate(r2);
       tmp32 = is_degenerate(s2);
       tmp32 = is_degenerate(sp2);
       tmp32 = is_degenerate(t2);
       tmp32 = is_degenerate(th2);
  (void) tmp32;


  typename R::Has_on_3 has_on
        = rep.has_on_3_object();
  bool tmp33a = has_on(l2,p2);
  bool tmp33b = has_on(t2,p2);
  bool tmp33c = has_on(tmp8,p2);
  (void) tmp33a;
  (void) tmp33b;
  (void) tmp33c;


  typename R::Has_on_bounded_side_3 has_on_bounded_side
        = rep.has_on_bounded_side_3_object();
  bool tmp34 = has_on_bounded_side(th2,p2);
  bool tmp34a = has_on_bounded_side(sp2,p2);
  bool tmp34b = has_on_bounded_side(iso1,p2);
  (void) tmp34;
  (void) tmp34a;
  (void) tmp34b;


  typename R::Has_on_unbounded_side_3 has_on_unbounded_side
        = rep.has_on_unbounded_side_3_object();
  bool tmp35 = has_on_unbounded_side(th2,p2);
  bool tmp35a = has_on_unbounded_side(sp2,p2);
  bool tmp35b = has_on_unbounded_side(iso1,p2);
  (void) tmp35;
  (void) tmp35a;
  (void) tmp35b;


  typename R::Has_on_boundary_3 has_on_boundary
        = rep.has_on_boundary_3_object();
  bool tmp36b = has_on_boundary(th2,p2);
  bool tmp36c = has_on_boundary(sp2,p2);
  bool tmp36d = has_on_boundary(iso1,p2);
  (void) tmp36b;
  (void) tmp36c;
  (void) tmp36d;


  typename R::Has_on_positive_side_3 has_on_positive_side
        = rep.has_on_positive_side_3_object();
  bool tmp37 = has_on_positive_side(h2,p2);
  (void) tmp37;


  typename R::Has_on_negative_side_3 has_on_negative_side
        = rep.has_on_negative_side_3_object();
  bool tmp38 = has_on_negative_side(h2,p2);
  (void) tmp38;


  typename R::Oriented_side_3 oriented_side
        = rep.oriented_side_3_object();
  Oriented_side tmp39 = oriented_side(h2,p2);
                tmp39 = oriented_side(sp9,p2);
  (void) tmp39;

  typename R::Bounded_side_3 bounded_side
        = rep.bounded_side_3_object();
  Bounded_side tmp39a = bounded_side(sp1,p2);
               tmp39a = bounded_side(th2,p2);
               tmp39a = bounded_side(iso1,p2);
  (void) tmp39a;

  typename R::Are_ordered_along_line_3 are_ordered_along_line
        = rep.are_ordered_along_line_3_object();
  bool tmp40 = are_ordered_along_line(p2,p3,p4);
  (void) tmp40;


  typename R::Are_strictly_ordered_along_line_3 are_strictly_ordered_along_line
        = rep.are_strictly_ordered_along_line_3_object();
  bool tmp41 = are_strictly_ordered_along_line(p2,p3,p4);
  (void) tmp41;


  typename R::Collinear_are_ordered_along_line_3 collinear_are_ordered_along_line
        = rep.collinear_are_ordered_along_line_3_object();
  bool tmp42 = collinear_are_ordered_along_line(p2,p2,p3);
  (void) tmp42;


  typename R::Collinear_are_strictly_ordered_along_line_3 collinear_are_strictly_ordered_along_line
        = rep.collinear_are_strictly_ordered_along_line_3_object();
  bool tmp43 = collinear_are_strictly_ordered_along_line(p2,p2,p3);
  (void) tmp43;


  typename R::Side_of_oriented_sphere_3 side_of_oriented_sphere
        = rep.side_of_oriented_sphere_3_object();
  Oriented_side tmp44 = side_of_oriented_sphere(p2,p3,p4,p5,p6);
  (void) tmp44;


  typename R::Side_of_bounded_sphere_3 side_of_bounded_sphere
        = rep.side_of_bounded_sphere_3_object();
  Bounded_side tmp45 = side_of_bounded_sphere(p2,p3,p4,p5,p6);
               tmp45 = side_of_bounded_sphere(p2,p3,p4,p6);
               tmp45 = side_of_bounded_sphere(p2,p3,p6);
  (void) tmp45;

  typename R::Angle_3 angle
        = rep.angle_3_object();
  Angle tmp46 = angle(p2,p3,p4);
  (void) tmp46;
  
  typename R::Construct_bisector_3 construct_bisector
        = rep.construct_bisector_3_object();
  Plane_3 tmp47 = construct_bisector(p2,p3);
  (void) tmp47;

  // kill warnings ...
  use(v1); use(d1); use(s1); use(l1); use(l4); use(l5);
  use(l6); use(d4); use(d5); use(d6); use(d7); use(r1);
  use(h1); use(h4); use(h5); use(h6); 
  use(sp4); use(sp5); use(sp6); use(sp7); use(sp8);
  use(t1); use(th1);
  use(tmp3); use(tmp3a);
  use(tmp9); use(tmp14a); use(tmp5); use(tmp6);
  use(tmp7); use(tmp71); use(sp1a); use(tmp72);
  use(tmp12a); use(tmp12aa); use(tmp12b);
  use(bb1); use(bb2); use(bb3); use(bb4); use(bb5); use(bb6); 
  use(r4); use(r5); use(l7); use(l8); use(v7); use(v8); use(v9); use(h8);
  use(cccit);

  // More tests, that require sqrt().
  typedef ::CGAL::Algebraic_structure_traits<FT> AST; 
  static const bool has_sqrt = 
      ! ::boost::is_same< ::CGAL::Null_functor, typename AST::Sqrt >::value;
  _test_new_3_sqrt(rep, ::CGAL::Boolean_tag<has_sqrt>() );
  
  return true;
}

#endif // CGAL_TEST_NEW_3_H

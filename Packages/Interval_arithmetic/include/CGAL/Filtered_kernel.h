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
// release       :
// release_date  :
//
// file          : include/CGAL/Filtered_kernel.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_FILTERED_KERNEL_H
#define CGAL_FILTERED_KERNEL_H

// This file contains the definition of a generic kernel filter.
//
// TODO:
// - at the moment, it's restricted to IA filtering, but this should be
//   generalized to allow static filters...
// - at the moment, only the predicates are filtered.
//   Constructions will come later.
// - proper treatment of constructive predicates remains a question.
// - the kernel only works with traits only and as a pure traits only.
// - split in several files.

#include <CGAL/basic.h>
#include <CGAL/Filter_predicate.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>

CGAL_BEGIN_NAMESPACE

// CK = construction kernel.
// EK = exact kernel called when needed by the filter.
// FK = filtering kernel
template <class CK,
          class EK,
          class FK = Simple_cartesian<Interval_nt_advanced>,
	  class C2E_Converter = Cartesian_converter<CK, EK>,
	  class C2F_Converter = Cartesian_converter<CK, FK,
	                Interval_converter<CGAL_TYPENAME_MSVC_NULL CK::RT> > >
class Filtered_kernel
{
public:
    // What to do with the tag ?
    // Probably this should not exist, should it ?
    // struct filter_tag{};
    // typedef filter_tag                                     Kernel_tag;
    // typedef typename CK::Kernel_tag                       Kernel_tag;
    typedef typename CK::Rep_tag                          Rep_tag;
    typedef typename CK::RT                               RT;
    typedef typename CK::FT                               FT;

    // Macros to define the types, predicates and constructions.

    // If we adopt this solution, then we can simply derive ?
#define CGAL_Filter_type(X) \
    typedef typename CK::X X##_base; \
    typedef typename CK::X X;

#define CGAL_Filter_pred(P, Pf) \
    typedef Filter_predicate<typename EK::P, typename FK::P, \
	                     C2E_Converter, C2F_Converter> P; \
    P Pf() const { return P(); }

#define CGAL_Filter_pred_type_only(P) \
    typedef Filter_predicate<typename EK::P, typename FK::P, \
	                     C2E_Converter, C2F_Converter> P;

    // The following could be used instead for some Cartesian predicates
    // that are exact : compare* and exact*.
#define CGAL_Filter_already_exact_pred(P, Pf) \
    typedef typename CK::P P; \
    P Pf() const { return P(); }

#define CGAL_Filter_cons(C, Cf) \
    typedef typename CK::C C; \
    C Cf() const { return C(); }

    // Types :

    // CGAL_Filter_type(RT) // ?
    // CGAL_Filter_type(FT) // ?

    CGAL_Filter_type(Point_2)
    CGAL_Filter_type(Vector_2)
    CGAL_Filter_type(Direction_2)
    CGAL_Filter_type(Segment_2)
    CGAL_Filter_type(Line_2)
    CGAL_Filter_type(Ray_2)
    CGAL_Filter_type(Triangle_2)
    CGAL_Filter_type(Circle_2)
    CGAL_Filter_type(Iso_rectangle_2)
    CGAL_Filter_type(Aff_transformation_2)
    // CGAL_Filter_type(Data_accessor_2) // ?
    // CGAL_Filter_type(Conic_2)         // ?

    CGAL_Filter_type(Point_3)
    CGAL_Filter_type(Vector_3)
    CGAL_Filter_type(Direction_3)
    CGAL_Filter_type(Segment_3)
    CGAL_Filter_type(Line_3)
    CGAL_Filter_type(Plane_3)
    CGAL_Filter_type(Ray_3)
    CGAL_Filter_type(Triangle_3)
    CGAL_Filter_type(Tetrahedron_3)
    CGAL_Filter_type(Sphere_3)
    CGAL_Filter_type(Iso_cuboid_3)
    CGAL_Filter_type(Aff_transformation_3)

    // Predicates and constructions :

    CGAL_Filter_cons(Construct_point_2, construct_point_2_object)
    CGAL_Filter_cons(Construct_vector_2, construct_vector_2_object)
    CGAL_Filter_cons(Construct_direction_2, construct_direction_2_object)
    CGAL_Filter_cons(Construct_segment_2, construct_segment_2_object)
    CGAL_Filter_cons(Construct_line_2, construct_line_2_object)
    CGAL_Filter_cons(Construct_ray_2, construct_ray_2_object)
    CGAL_Filter_cons(Construct_circle_2, construct_circle_2_object)
    CGAL_Filter_cons(Construct_triangle_2, construct_triangle_2_object)
    CGAL_Filter_cons(Construct_iso_rectangle_2,
	    construct_iso_rectangle_2_object)
    CGAL_Filter_cons(Construct_aff_transformation_2,
	    construct_aff_transformation_2_object)
    CGAL_Filter_cons(Construct_point_on_2, construct_point_on_2_object)
    CGAL_Filter_cons(Construct_second_point_on_2,
	    construct_second_point_on_2_object)
    CGAL_Filter_cons(Construct_source_point_2, construct_source_point_2_object)
    CGAL_Filter_cons(Construct_target_point_2, construct_target_point_2_object)
    CGAL_Filter_cons(Construct_min_point_2, construct_min_point_2_object)
    CGAL_Filter_cons(Construct_max_point_2, construct_max_point_2_object)
    CGAL_Filter_cons(Construct_direction_of_line_2,
	construct_direction_of_line_2_object)
    CGAL_Filter_cons(Construct_direction_of_ray_2,
	construct_direction_of_ray_2_object)
    CGAL_Filter_cons(Construct_supporting_line_2,
	construct_supporting_line_2_object)
    CGAL_Filter_cons(Construct_perpendicular_vector_2,
	construct_perpendicular_vector_2_object)
    CGAL_Filter_cons(Construct_perpendicular_direction_2,
	construct_perpendicular_direction_2_object)
    CGAL_Filter_cons(Construct_perpendicular_line_2,
	construct_perpendicular_line_2_object)
    CGAL_Filter_cons(Construct_midpoint_2, construct_midpoint_2_object)
    CGAL_Filter_cons(Construct_center_2, construct_center_2_object)
    CGAL_Filter_cons(Construct_circumcenter_2, construct_circumcenter_2_object)
    CGAL_Filter_cons(Construct_bisector_2, construct_bisector_2_object)
    CGAL_Filter_cons(Construct_opposite_segment_2,
	construct_opposite_segment_2_object)
    CGAL_Filter_cons(Construct_opposite_ray_2, construct_opposite_ray_2_object)
    CGAL_Filter_cons(Construct_opposite_line_2,
	    construct_opposite_line_2_object)
    CGAL_Filter_cons(Construct_opposite_triangle_2,
	construct_opposite_triangle_2_object)
    CGAL_Filter_cons(Construct_opposite_circle_2,
	    construct_opposite_circle_2_object)
    CGAL_Filter_cons(Assign_2, assign_2_object)
    CGAL_Filter_cons(Transform_2, transform_2_object)
    CGAL_Filter_cons(Intersect_2, intersect_2_object)
    CGAL_Filter_cons(Compute_y_at_x_2, compute_y_at_x_2_object)
    CGAL_Filter_cons(Compute_squared_distance_2,
	    Compute_squared_distance_2_object)
    CGAL_Filter_cons(Compute_squared_length_2,
	    Compute_squared_length_2_object)
    CGAL_Filter_cons(Compute_squared_radius_2, Compute_squared_radius_2_object)


    CGAL_Filter_pred(Equal_2, equal_2_object)
    CGAL_Filter_pred(Equal_x_2, equal_x_2_object)
    CGAL_Filter_pred(Equal_y_2, equal_y_2_object)
    CGAL_Filter_pred(Equal_xy_2, equal_xy_2_object)
    CGAL_Filter_pred(Less_x_2, less_x_2_object)
    CGAL_Filter_pred(Less_y_2, less_y_2_object)
    CGAL_Filter_pred(Less_xy_2, less_xy_2_object)
    CGAL_Filter_pred(Less_yx_2, less_yx_2_object)
    CGAL_Filter_pred(Compare_x_2, compare_x_2_object)
    CGAL_Filter_pred(Compare_y_2, compare_y_2_object)
    CGAL_Filter_pred(Compare_xy_2, compare_xy_2_object)
    CGAL_Filter_pred(Compare_y_at_x_2, compare_y_at_x_2_object)
    CGAL_Filter_pred(Counterclockwise_in_between_2,
	counterclockwise_in_between_2_object)
    CGAL_Filter_pred(Leftturn_2, leftturn_2_object)
    CGAL_Filter_pred(Collinear_2, collinear_2_object)
    CGAL_Filter_pred(Orientation_2, orientation_2_object)
    CGAL_Filter_pred(Side_of_oriented_circle_2,
	    side_of_oriented_circle_2_object)
    CGAL_Filter_pred(Side_of_bounded_circle_2, side_of_bounded_circle_2_object)
    CGAL_Filter_pred(Is_horizontal_2, is_horizontal_2_object)
    CGAL_Filter_pred(Is_vertical_2, is_vertical_2_object)
    CGAL_Filter_pred(Is_degenerate_2, is_degenerate_2_object)
    CGAL_Filter_pred(Has_on_2, has_on_2_object)
    CGAL_Filter_pred(Collinear_has_on_2, collinear_has_on_2_object)
    CGAL_Filter_pred(Has_on_bounded_side_2, has_on_bounded_side_2_object)
    CGAL_Filter_pred(Has_on_unbounded_side_2, has_on_unbounded_side_2_object)
    CGAL_Filter_pred(Has_on_boundary_2, has_on_boundary_2_object)
    CGAL_Filter_pred(Has_on_positive_side_2, has_on_positive_side_2_object)
    CGAL_Filter_pred(Has_on_negative_side_2, has_on_negative_side_2_object)
    CGAL_Filter_pred(Oriented_side_2, oriented_side_2_object)
    CGAL_Filter_pred(Are_ordered_along_line_2, are_ordered_along_line_2_object)
    CGAL_Filter_pred(Are_strictly_ordered_along_line_2,
	are_strictly_ordered_along_line_2_object)
    CGAL_Filter_pred(Collinear_are_ordered_along_line_2,
	collinear_are_ordered_along_line_2_object)
    CGAL_Filter_pred(Collinear_are_strictly_ordered_along_line_2,
	collinear_are_strictly_ordered_along_line_2_object)


    CGAL_Filter_cons(Construct_point_3, construct_point_3_object)
    CGAL_Filter_cons(Construct_vector_3, construct_vector_3_object)
    CGAL_Filter_cons(Construct_direction_3, construct_direction_3_object)
    CGAL_Filter_cons(Construct_segment_3, construct_segment_3_object)
    CGAL_Filter_cons(Construct_plane_3, construct_plane_3_object)
    CGAL_Filter_cons(Construct_line_3, construct_line_3_object)
    CGAL_Filter_cons(Construct_ray_3, construct_ray_3_object)
    CGAL_Filter_cons(Construct_sphere_3, construct_sphere_3_object)
    CGAL_Filter_cons(Construct_triangle_3, construct_triangle_3_object)
    CGAL_Filter_cons(Construct_tetrahedron_3, construct_tetrahedron_3_object)
    CGAL_Filter_cons(Construct_iso_cuboid_3, construct_iso_cuboid_3_object)
    CGAL_Filter_cons(Construct_aff_transformation_3,
	construct_aff_transformation_3_object)
    CGAL_Filter_cons(Construct_point_on_3, construct_point_on_3_object)
    CGAL_Filter_cons(Construct_second_point_on_3,
	    construct_second_point_on_3_object)
    CGAL_Filter_cons(Construct_source_point_3, construct_source_point_3_object)
    CGAL_Filter_cons(Construct_target_point_3, construct_target_point_3_object)
    CGAL_Filter_cons(Construct_min_point_3, construct_min_point_3_object)
    CGAL_Filter_cons(Construct_max_point_3, construct_max_point_3_object)
    CGAL_Filter_cons(Construct_direction_of_line_3, 
	construct_direction_of_line_3_object)
    CGAL_Filter_cons(Construct_direction_of_ray_3,
	construct_direction_of_ray_3_object)
    CGAL_Filter_cons(Construct_supporting_line_3,
	    construct_supporting_line_3_object)
    CGAL_Filter_cons(Construct_perpendicular_plane_3,
	construct_perpendicular_plane_3_object)
    CGAL_Filter_cons(Construct_perpendicular_line_3,
	construct_perpendicular_line_3_object)
    CGAL_Filter_cons(Construct_midpoint_3, construct_midpoint_3_object)
    CGAL_Filter_cons(Construct_center_3, construct_center_3_object)
    CGAL_Filter_cons(Construct_circumcenter_3, construct_circumcenter_3_object)
    CGAL_Filter_cons(Construct_cross_product_vector_3,
	    construct_cross_product_vector_3_object)
    CGAL_Filter_cons(Construct_opposite_segment_3,
	construct_opposite_segment_3_object)
    CGAL_Filter_cons(Construct_opposite_ray_3, construct_opposite_ray_3_object)
    CGAL_Filter_cons(Construct_opposite_line_3,
	    construct_opposite_line_3_object)
    CGAL_Filter_cons(Construct_opposite_plane_3,
	    construct_opposite_plane_3_object)
    CGAL_Filter_cons(Construct_supporting_plane_3,
	construct_supporting_plane_3_object)
    CGAL_Filter_cons(Transform_3, transform_3_object)
    CGAL_Filter_cons(Assign_3, assign_3_object)
    CGAL_Filter_cons(Intersect_3, intersect_3_object)
    CGAL_Filter_cons(Compute_squared_distance_3,
	    compute_squared_distance_3_object)
    CGAL_Filter_cons(Compute_squared_length_3, compute_squared_length_3_object)
    CGAL_Filter_cons(Compute_squared_radius_3, compute_squared_radius_3_object)


    CGAL_Filter_pred(Equal_3, equal_3_object)
    CGAL_Filter_pred(Equal_x_3, equal_x_3_object)
    CGAL_Filter_pred(Equal_y_3, equal_y_3_object)
    CGAL_Filter_pred(Equal_z_3, equal_z_3_object)
    CGAL_Filter_pred(Equal_xy_3, equal_xy_3_object)
    CGAL_Filter_pred(Equal_xyz_3, equal_xyz_3_object)
    CGAL_Filter_pred(Less_x_3, less_x_3_object)
    CGAL_Filter_pred(Less_y_3, less_y_3_object)
    CGAL_Filter_pred(Less_z_3, less_z_3_object)
    CGAL_Filter_pred(Less_xy_3, less_xy_3_object)
    CGAL_Filter_pred(Less_xyz_3, less_xyz_3_object)
    CGAL_Filter_pred(Compare_x_3, compare_x_3_object)
    CGAL_Filter_pred(Compare_y_3, compare_y_3_object)
    CGAL_Filter_pred(Compare_z_3, compare_z_3_object)
    CGAL_Filter_pred(Compare_xy_3, compare_xy_3_object)
    CGAL_Filter_pred(Compare_xyz_3, compare_xyz_3_object)
    CGAL_Filter_pred(Collinear_3, collinear_3_object)
    CGAL_Filter_pred(Coplanar_3, coplanar_3_object)
    CGAL_Filter_pred(Coplanar_orientation_3, coplanar_orientation_3_object)
    CGAL_Filter_pred(Coplanar_side_of_bounded_circle_3,
	    coplanar_side_of_bounded_circle_3_object)
    CGAL_Filter_pred(Orientation_3, orientation_3_object)
    CGAL_Filter_pred(Is_degenerate_3, is_degenerate_3_object)
    CGAL_Filter_pred(Has_on_3, has_on_3_object)
    CGAL_Filter_pred(Has_on_bounded_side_3, has_on_bounded_side_3_object)
    CGAL_Filter_pred(Has_on_unbounded_side_3, has_on_unbounded_side_3_object)
    CGAL_Filter_pred(Has_on_boundary_3, has_on_boundary_3_object)
    CGAL_Filter_pred(Has_on_positive_side_3, has_on_positive_side_3_object)
    CGAL_Filter_pred(Has_on_negative_side_3, has_on_negative_side_3_object)
    CGAL_Filter_pred(Oriented_side_3, oriented_side_3_object)
    CGAL_Filter_pred(Are_ordered_along_line_3, are_ordered_along_line_3_object)
    CGAL_Filter_pred(Are_strictly_ordered_along_line_3,
	    are_strictly_ordered_along_line_3_object)
    CGAL_Filter_pred(Collinear_are_ordered_along_line_3,
	    collinear_are_ordered_along_line_3_object)
    CGAL_Filter_pred(Collinear_are_strictly_ordered_along_line_3,
	    collinear_are_strictly_ordered_along_line_3)
    CGAL_Filter_pred(Side_of_oriented_sphere_3,
	    side_of_oriented_sphere_3_object)
    CGAL_Filter_pred(Side_of_bounded_sphere_3, side_of_bounded_sphere_3_object)

    // Constructive predicates
    // They are problematic, as their constructor will systematically
    // construct the _exact_ object, and thus will be uselessly costly.

    CGAL_Filter_pred_type_only(Compare_distance_to_point_2)
    CGAL_Filter_pred_type_only(Less_distance_to_point_2)
    CGAL_Filter_pred_type_only(Less_signed_distance_to_line_2)
    CGAL_Filter_pred_type_only(Less_rotate_ccw_2)
    CGAL_Filter_pred_type_only(Left_of_line_2)
    CGAL_Filter_pred_type_only(Compare_distance_to_point_3)
    CGAL_Filter_pred_type_only(Less_distance_to_point_3)

    Less_distance_to_point_2
    less_distance_to_point_2_object(const Point_2& p) const
    { return Less_distance_to_point_2(p); }

    Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object(const Point_2& p,
                                          const Point_2& q) const
    { return Less_signed_distance_to_line_2(p,q); }

    Less_rotate_ccw_2
    less_rotate_ccw_2_object(const Point_2& p) const
    { return Less_rotate_ccw_2(p); }

    Left_of_line_2
    left_of_line_2_object(const Point_2& p, const Point_2& q) const
    { return Left_of_line_2(p,q); }

    Less_distance_to_point_3
    less_distance_to_point_3_object(const Point_3& p) const
    { return Less_distance_to_point_3(p); }

    Compare_distance_to_point_3
    compare_distance_to_point_3_object(const Point_3& p) const
    { return Compare_distance_to_point_3(p); }

    // CGAL_Filter_cons(Construct_point_d, construct_point_d_object)
};

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_KERNEL_H

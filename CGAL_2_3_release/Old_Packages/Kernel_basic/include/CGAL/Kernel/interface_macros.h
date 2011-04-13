// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
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
// file          : include/CGAL/Kernel/interface_macros.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Sylvain Pion, Susan Hert
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred and CGAL_Kernel_cons.
// And they are #undefed at the end of this file.
// Note : there are also CGAL_Kernel_pred2 and CGAL_Kernel_cons2 for those
// which first argument contains a comma.

CGAL_Kernel_cons(CGALi::Construct<Point_2>,
	         Construct_point_2,
		 construct_point_2_object)
CGAL_Kernel_cons(CGALi::Construct<Vector_2>,
	         Construct_vector_2,
		 construct_vector_2_object)
CGAL_Kernel_cons(CGALi::Construct<Direction_2>,
	         Construct_direction_2,
		 construct_direction_2_object)
CGAL_Kernel_cons(CGALi::Construct<Segment_2>,
	         Construct_segment_2,
		 construct_segment_2_object)
CGAL_Kernel_cons(CGALi::Construct<Line_2>,
	         Construct_line_2,
		 construct_line_2_object)
CGAL_Kernel_cons(CGALi::Construct<Ray_2>,
	         Construct_ray_2,
		 construct_ray_2_object)
CGAL_Kernel_cons(CGALi::Construct<Circle_2>,
	         Construct_circle_2,
		 construct_circle_2_object)
CGAL_Kernel_cons(CGALi::Construct<Triangle_2>,
	         Construct_triangle_2,
		 construct_triangle_2_object)
CGAL_Kernel_cons(CGALi::Construct<Iso_rectangle_2>,
	         Construct_iso_rectangle_2,
		 construct_iso_rectangle_2_object)
CGAL_Kernel_cons(CGALi::Call_make_object_to_get<Object_2>,
                 Construct_object_2,
                 construct_object_2_object)

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Construct<Aff_transformation_2>,
	         Construct_aff_transformation_2,
		 construct_aff_transformation_2_object)
#endif // CGAL_NO_DEPRECATED_CODE

CGAL_Kernel_cons(CGALi::Call_point_to_get<Point_2>,
	         Construct_point_on_2,
		 construct_point_on_2_object)

CGAL_Kernel_cons(CGALi::Call_projection_to_get<Point_2>,
	         Construct_projected_point_2,
		 construct_projected_point_2_object)

CGAL_Kernel_cons(CGALi::Construct_projected_xy_point<Point_2>,
	         Construct_projected_xy_point_2,
		 construct_projected_xy_point_2_object)

CGAL_Kernel_cons(CGALi::Construct_scaled_vector<Vector_2>,
	         Construct_scaled_vector_2,
		 construct_scaled_vector_2_object)

CGAL_Kernel_cons(CGALi::Construct_translated_point<Point_2>,
	         Construct_translated_point_2,
		 construct_translated_point_2_object)

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_second_point_to_get<Point_2>,
	         Construct_second_point_on_2,
		 construct_second_point_on_2_object)
CGAL_Kernel_cons(CGALi::Call_source_to_get<Point_2>,
	         Construct_source_point_2,
		 construct_source_point_2_object)
CGAL_Kernel_cons(CGALi::Call_target_to_get<Point_2>,
	         Construct_target_point_2,
		 construct_target_point_2_object)
CGAL_Kernel_cons(CGALi::Call_min_to_get<Point_2>,
	         Construct_min_point_2,
		 construct_min_point_2_object)
CGAL_Kernel_cons(CGALi::Call_max_to_get<Point_2>,
	         Construct_max_point_2,
		 construct_max_point_2_object)
#endif // CGAL_NO_DEPRECATED_CODE

CGAL_Kernel_cons(CGALi::Call_vertex_to_get<Point_2>,
	         Construct_vertex_2,
		 construct_vertex_2_object)

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_direction_to_get<Direction_2>,
	         Construct_direction_of_line_2,
		 construct_direction_of_line_2_object)
CGAL_Kernel_cons(CGALi::Call_direction_to_get<Direction_2>,
	         Construct_direction_of_ray_2,
		 construct_direction_of_ray_2_object)
#endif // CGAL_NO_DEPRECATED_CODE

CGAL_Kernel_cons(CGALi::Call_supporting_line_to_get<Line_2>,
	         Construct_supporting_line_2,
		 construct_supporting_line_2_object)
CGAL_Kernel_cons(CGALi::Call_perpendicular_to_get<Vector_2>,
	         Construct_perpendicular_vector_2,
		 construct_perpendicular_vector_2_object)
CGAL_Kernel_cons(CGALi::Call_perpendicular_to_get<Direction_2>,
	         Construct_perpendicular_direction_2,
		 construct_perpendicular_direction_2_object)
CGAL_Kernel_cons(CGALi::Call_perpendicular_to_get<Line_2>,
	         Construct_perpendicular_line_2,
		 construct_perpendicular_line_2_object)
CGAL_Kernel_cons(CGALi::p_Midpoint<Point_2>,
	         Construct_midpoint_2,
		 construct_midpoint_2_object)
CGAL_Kernel_cons(CGALi::p_Center<Point_2>,
	         Construct_center_2,
		 construct_center_2_object)
CGAL_Kernel_cons(CGALi::p_Circumcenter<Point_2>,
	         Construct_circumcenter_2,
		 construct_circumcenter_2_object)
CGAL_Kernel_cons(CGALi::p_Centroid<Point_2>,
	         Construct_centroid_2,
		 construct_centroid_2_object)
CGAL_Kernel_cons2(CGALi::pl_Bisector<Point_2, Line_2>,
	         Construct_bisector_2,
		 construct_bisector_2_object)
CGAL_Kernel_cons(CGALi::v_Opposite<Direction_2>,
	         Construct_opposite_direction_2,
		 construct_opposite_direction_2_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Segment_2>,
	         Construct_opposite_segment_2,
		 construct_opposite_segment_2_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Ray_2>,
	         Construct_opposite_ray_2,
		 construct_opposite_ray_2_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Line_2>,
	         Construct_opposite_line_2,
		 construct_opposite_line_2_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Triangle_2>,
	         Construct_opposite_triangle_2,
		 construct_opposite_triangle_2_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Circle_2>,
	         Construct_opposite_circle_2,
		 construct_opposite_circle_2_object)
CGAL_Kernel_cons(CGALi::v_Opposite<Vector_2>,
	         Construct_opposite_vector_2,
		 construct_opposite_vector_2_object)
CGAL_Kernel_cons(CGALi::Assign,
	         Assign_2,
		 assign_2_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_transform,
	         Transform_2,
		 transform_2_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Intersect,
	         Intersect_2,
		 intersect_2_object)
CGAL_Kernel_pred(CGALi::Do_intersect,
	         Do_intersect_2,
		 do_intersect_2_object)

CGAL_Kernel_cons(CGALi::Call_y_at_x_to_get<FT>,
	         Compute_y_at_x_2,
		 compute_y_at_x_2_object)
CGAL_Kernel_cons(CGALi::Call_squared_distance<FT>,
	         Compute_squared_distance_2,
		 compute_squared_distance_2_object)
CGAL_Kernel_cons(CGALi::Call_squared_length_to_get<FT>,
	         Compute_squared_length_2,
		 compute_squared_length_2_object)
CGAL_Kernel_cons(CGALi::Call_squared_radius<FT>,
	         Compute_squared_radius_2,
		 compute_squared_radius_2_object)
CGAL_Kernel_cons(CGALi::Call_area_to_get<FT>,
	         Compute_area_2,
		 compute_area_2_object)
CGAL_Kernel_pred(CGALi::Equal,
	         Equal_2,
		 equal_2_object)
CGAL_Kernel_pred(CGALi::Equal_x,
	         Equal_x_2,
		 equal_x_2_object)
CGAL_Kernel_pred(CGALi::Equal_y,
	         Equal_y_2,
		 equal_y_2_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_pred(CGALi::Equal_xy,
	         Equal_xy_2,
		 equal_xy_2_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_pred(CGALi::Less_x,
	         Less_x_2,
		 less_x_2_object)
CGAL_Kernel_pred(CGALi::Less_y,
	         Less_y_2,
		 less_y_2_object)
CGAL_Kernel_pred(CGAL::p_Less_xy<Point_2>,
	         Less_xy_2,
		 less_xy_2_object)
CGAL_Kernel_pred(CGAL::p_Less_yx<Point_2>,
	         Less_yx_2,
		 less_yx_2_object)
CGAL_Kernel_pred(CGALi::Compare_x,
	         Compare_x_2,
		 compare_x_2_object)
CGAL_Kernel_pred(CGALi::Compare_y,
	         Compare_y_2,
		 compare_y_2_object)
CGAL_Kernel_pred(CGALi::Compare_xy,
	         Compare_xy_2,
		 compare_xy_2_object)
CGAL_Kernel_pred(CGALi::Compare_y_at_x,
	         Compare_y_at_x_2,
		 compare_y_at_x_2_object)
CGAL_Kernel_pred(CGALi::Compare_x_at_y,
	         Compare_x_at_y_2,
		 compare_x_at_y_2_object)
CGAL_Kernel_pred(CGALi::Compare_distance<Point_2>,
	         Compare_distance_2,
		 compare_distance_2_object)
CGAL_Kernel_pred(CGAL ::p_Less_dist_to_point<Point_2>,
                 Less_distance_to_point_2,
                 less_distance_to_point_2_object)
CGAL_Kernel_pred(CGAL ::p_Less_dist_to_line_2<Point_2>,
                 Less_signed_distance_to_line_2,
                 less_signed_distance_to_line_2_object)
CGAL_Kernel_pred(CGAL ::p_Less_rotate_ccw<Point_2>,
	         Less_rotate_ccw_2,
		 less_rotate_ccw_2_object)
CGAL_Kernel_pred(CGALi::Compare_angle_with_x_axis<Direction_2>,
	         Compare_angle_with_x_axis_2,
		 compare_angle_with_x_axis_2_object)
CGAL_Kernel_pred(CGALi::Counterclockwise_in_between,
	         Counterclockwise_in_between_2,
		 counterclockwise_in_between_2_object)
CGAL_Kernel_pred(CGAL ::p_Left_turn<Point_2>,
                 Left_turn_2,
                 left_turn_2_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_pred(CGAL ::p_Left_turn<Point_2>,
	         Leftturn_2,
		 leftturn_2_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_pred(CGALi::p_Angle,
	         Angle_2,
		 angle_2_object)
CGAL_Kernel_pred(CGALi::Collinear,
	         Collinear_2,
		 collinear_2_object)
CGAL_Kernel_pred(CGAL ::p_Orientation<Point_2>,
	         Orientation_2,
		 orientation_2_object)
CGAL_Kernel_pred(CGALi::Side_of_oriented_circle,
	         Side_of_oriented_circle_2,
		 side_of_oriented_circle_2_object)
CGAL_Kernel_pred(CGALi::Side_of_bounded_circle,
	         Side_of_bounded_circle_2,
		 side_of_bounded_circle_2_object)
CGAL_Kernel_pred(CGALi::Call_is_horizontal,
	         Is_horizontal_2,
		 is_horizontal_2_object)
CGAL_Kernel_pred(CGALi::Call_is_vertical,
	         Is_vertical_2,
		 is_vertical_2_object)
CGAL_Kernel_pred(CGALi::Call_is_degenerate,
	         Is_degenerate_2,
		 is_degenerate_2_object)
CGAL_Kernel_pred(CGALi::Call_has_on,
	         Has_on_2,
		 has_on_2_object)
CGAL_Kernel_pred(CGALi::Call_collinear_has_on,
	         Collinear_has_on_2,
		 collinear_has_on_2_object)
CGAL_Kernel_pred(CGALi::Call_has_on_bounded_side,
	         Has_on_bounded_side_2,
		 has_on_bounded_side_2_object)
CGAL_Kernel_pred(CGALi::Call_has_on_unbounded_side,
	         Has_on_unbounded_side_2,
		 has_on_unbounded_side_2_object)
CGAL_Kernel_pred(CGALi::Call_has_on_boundary,
	         Has_on_boundary_2,
		 has_on_boundary_2_object)
CGAL_Kernel_pred(CGALi::Call_has_on_positive_side,
	         Has_on_positive_side_2,
		 has_on_positive_side_2_object)
CGAL_Kernel_pred(CGALi::Call_has_on_negative_side,
	         Has_on_negative_side_2,
		 has_on_negative_side_2_object)
CGAL_Kernel_pred(CGALi::Call_oriented_side,
	         Oriented_side_2,
		 oriented_side_2_object)
CGAL_Kernel_pred(CGALi::Call_bounded_side,
	         Bounded_side_2,
		 bounded_side_2_object)
CGAL_Kernel_pred(CGALi::Are_ordered_along_line,
	         Are_ordered_along_line_2,
		 are_ordered_along_line_2_object)
CGAL_Kernel_pred(CGALi::Are_strictly_ordered_along_line,
	         Are_strictly_ordered_along_line_2,
		 are_strictly_ordered_along_line_2_object)
CGAL_Kernel_pred(CGALi::Collinear_are_ordered_along_line,
	         Collinear_are_ordered_along_line_2,
		 collinear_are_ordered_along_line_2_object)
CGAL_Kernel_pred(CGALi::Collinear_are_strictly_ordered_along_line,
	         Collinear_are_strictly_ordered_along_line_2,
		 collinear_are_strictly_ordered_along_line_2_object)

CGAL_Kernel_cons(CGALi::Construct<Point_3>,
	         Construct_point_3,
		 construct_point_3_object)
CGAL_Kernel_cons(CGALi::Construct<Vector_3>,
	         Construct_vector_3,
		 construct_vector_3_object)
CGAL_Kernel_cons(CGALi::Construct<Direction_3>,
	         Construct_direction_3,
		 construct_direction_3_object)
CGAL_Kernel_cons(CGALi::Construct<Segment_3>,
	         Construct_segment_3,
		 construct_segment_3_object)
CGAL_Kernel_cons(CGALi::Construct<Plane_3>,
	         Construct_plane_3,
		 construct_plane_3_object)
CGAL_Kernel_cons(CGALi::Construct<Line_3>,
	         Construct_line_3,
		 construct_line_3_object)
CGAL_Kernel_cons(CGALi::Construct<Ray_3>,
	         Construct_ray_3,
		 construct_ray_3_object)
CGAL_Kernel_cons(CGALi::Construct<Sphere_3>,
	         Construct_sphere_3,
		 construct_sphere_3_object)
CGAL_Kernel_cons(CGALi::Construct<Triangle_3>,
	         Construct_triangle_3,
		 construct_triangle_3_object)
CGAL_Kernel_cons(CGALi::Construct<Tetrahedron_3>,
	         Construct_tetrahedron_3,
		 construct_tetrahedron_3_object)
CGAL_Kernel_cons(CGALi::Construct<Iso_cuboid_3>,
	         Construct_iso_cuboid_3,
		 construct_iso_cuboid_3_object)
CGAL_Kernel_cons(CGALi::Call_make_object_to_get<Object_3>,
                 Construct_object_3,
                 construct_object_3_object)

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Construct<Aff_transformation_3>,
	         Construct_aff_transformation_3,
		 construct_aff_transformation_3_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_point_to_get<Point_3>,
	         Construct_point_on_3,
		 construct_point_on_3_object)
CGAL_Kernel_cons(CGALi::Call_projection_to_get<Point_3>,
	         Construct_projected_point_3,
		 construct_projected_point_3_object)
CGAL_Kernel_cons(CGALi::p_Lifted<Point_3>,
	         Construct_lifted_point_3,
		 construct_lifted_point_3_object)
CGAL_Kernel_cons(CGALi::Construct_scaled_vector<Vector_3>,
	         Construct_scaled_vector_3,
		 construct_scaled_vector_3_object)
CGAL_Kernel_cons(CGALi::Construct_translated_point<Point_3>,
	         Construct_translated_point_3,
		 construct_translated_point_3_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_second_point_to_get<Point_3>,
	         Construct_second_point_on_3,
		 construct_second_point_on_3_object)
CGAL_Kernel_cons(CGALi::Call_source_to_get<Point_3>,
	         Construct_source_point_3,
		 construct_source_point_3_object)
CGAL_Kernel_cons(CGALi::Call_target_to_get<Point_3>,
	         Construct_target_point_3,
		 construct_target_point_3_object)
CGAL_Kernel_cons(CGALi::Call_min_to_get<Point_3>,
	         Construct_min_point_3,
		 construct_min_point_3_object)
CGAL_Kernel_cons(CGALi::Call_max_to_get<Point_3>,
	         Construct_max_point_3,
		 construct_max_point_3_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_vertex_to_get<Point_3>,
	         Construct_vertex_3,
		 construct_vertex_3_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_direction_to_get<Direction_3>,
	         Construct_direction_of_line_3,
		 construct_direction_of_line_3_object)
CGAL_Kernel_cons(CGALi::Call_direction_to_get<Direction_3>,
	         Construct_direction_of_ray_3,
		 construct_direction_of_ray_3_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_supporting_line_to_get<Line_3>,
	         Construct_supporting_line_3,
		 construct_supporting_line_3_object)
CGAL_Kernel_cons(CGALi::Call_perpendicular_plane_to_get<Plane_3>,
	         Construct_perpendicular_plane_3,
		 construct_perpendicular_plane_3_object)
CGAL_Kernel_cons(CGALi::Call_perpendicular_line_to_get<Line_3>,
	         Construct_perpendicular_line_3,
		 construct_perpendicular_line_3_object)
CGAL_Kernel_cons(CGALi::Call_orthogonal_vector_to_get<Vector_3>,
	         Construct_orthogonal_vector_3,
		 construct_orthogonal_vector_3_object)
CGAL_Kernel_cons(CGALi::v_Base<Vector_3>,
	         Construct_base_vector_3,
		 construct_base_vector_3_object)
CGAL_Kernel_cons(CGALi::p_Midpoint<Point_3>,
	         Construct_midpoint_3,
		 construct_midpoint_3_object)
CGAL_Kernel_cons(CGALi::p_Center<Point_3>,
	         Construct_center_3,
		 construct_center_3_object)
CGAL_Kernel_cons(CGALi::p_Circumcenter<Point_3>,
	         Construct_circumcenter_3,
		 construct_circumcenter_3_object)
CGAL_Kernel_cons(CGALi::p_Centroid<Point_3>,
	         Construct_centroid_3,
		 construct_centroid_3_object)
CGAL_Kernel_cons(CGALi::v_Cross_product<Vector_3>,
	         Construct_cross_product_vector_3,
		 construct_cross_product_vector_3_object)
CGAL_Kernel_cons(CGALi::v_Opposite<Direction_3>,
	         Construct_opposite_direction_3,
		 construct_opposite_direction_3_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Segment_3>,
	         Construct_opposite_segment_3,
		 construct_opposite_segment_3_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Ray_3>,
	         Construct_opposite_ray_3,
		 construct_opposite_ray_3_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Line_3>,
	         Construct_opposite_line_3,
		 construct_opposite_line_3_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Plane_3>,
	         Construct_opposite_plane_3,
		 construct_opposite_plane_3_object)
CGAL_Kernel_cons(CGALi::Call_opposite_to_get<Sphere_3>,
	         Construct_opposite_sphere_3,
		 construct_opposite_sphere_3_object)
CGAL_Kernel_cons(CGALi::v_Opposite<Vector_3>,
	         Construct_opposite_vector_3,
		 construct_opposite_vector_3_object)
CGAL_Kernel_cons(CGALi::Call_supporting_plane_to_get<Plane_3>,
	         Construct_supporting_plane_3,
		 construct_supporting_plane_3_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Call_transform,
	         Transform_3,
		 transform_3_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_cons(CGALi::Assign,
	         Assign_3,
		 assign_3_object)
CGAL_Kernel_cons(CGALi::Intersect,
	         Intersect_3,
		 intersect_3_object)
CGAL_Kernel_cons(CGALi::Do_intersect,
	         Do_intersect_3,
		 do_intersect_3_object)
CGAL_Kernel_cons(CGALi::Call_squared_distance<FT>,
	         Compute_squared_distance_3,
		 compute_squared_distance_3_object)
CGAL_Kernel_cons(CGALi::Call_squared_length_to_get<FT>,
	         Compute_squared_length_3,
		 compute_squared_length_3_object)
CGAL_Kernel_cons(CGALi::Call_squared_radius<FT>,
	         Compute_squared_radius_3,
		 compute_squared_radius_3_object)
CGAL_Kernel_cons(CGALi::Call_squared_area_to_get<FT>,
	         Compute_squared_area_3,
		 compute_squared_area_3_object)
CGAL_Kernel_cons(CGALi::Call_volume_to_get<FT>,
	         Compute_volume_3,
		 compute_volume_3_object)
CGAL_Kernel_pred(CGALi::Equal,
	         Equal_3,
		 equal_3_object)
CGAL_Kernel_pred(CGALi::Equal_x,
	         Equal_x_3,
		 equal_x_3_object)
CGAL_Kernel_pred(CGALi::Equal_y,
	         Equal_y_3,
		 equal_y_3_object)
CGAL_Kernel_pred(CGALi::Equal_z,
	         Equal_z_3,
		 equal_z_3_object)
CGAL_Kernel_pred(CGALi::Equal_xy,
	         Equal_xy_3,
		 equal_xy_3_object)
#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_pred(CGALi::Equal_xyz,
	         Equal_xyz_3,
		 equal_xyz_3_object)
#endif // CGAL_NO_DEPRECATED_CODE
CGAL_Kernel_pred(CGALi::Less_x,
	         Less_x_3,
		 less_x_3_object)
CGAL_Kernel_pred(CGALi::Less_y,
	         Less_y_3,
		 less_y_3_object)
CGAL_Kernel_pred(CGALi::Less_z,
	         Less_z_3,
		 less_z_3_object)
CGAL_Kernel_pred(CGAL::p_Less_xy<Point_3>,
	         Less_xy_3,
		 less_xy_3_object)
CGAL_Kernel_pred(CGALi::Less_xyz,
	         Less_xyz_3,
		 less_xyz_3_object)
CGAL_Kernel_pred(CGALi::Compare_x,
	         Compare_x_3,
		 compare_x_3_object)
CGAL_Kernel_pred(CGALi::Compare_y,
	         Compare_y_3,
		 compare_y_3_object)
CGAL_Kernel_pred(CGALi::Compare_z,
	         Compare_z_3,
		 compare_z_3_object)
CGAL_Kernel_pred(CGALi::Compare_xy,
	         Compare_xy_3,
		 compare_xy_3_object)
CGAL_Kernel_pred(CGALi::Compare_xyz,
	         Compare_xyz_3,
		 compare_xyz_3_object)
CGAL_Kernel_pred(CGALi::Compare_distance<Point_3>,
	         Compare_distance_3,
		 compare_distance_3_object)
CGAL_Kernel_pred(CGAL ::p_Less_dist_to_point<Point_3>,
                 Less_distance_to_point_3,
                 less_distance_to_point_3_object)
CGAL_Kernel_pred2(CGALi::Less_signed_distance_to_plane<Plane_3, Point_3>,
	         Less_signed_distance_to_plane_3,
		 less_signed_distance_to_plane_3_object)
CGAL_Kernel_pred(CGALi::p_Angle,
	         Angle_3,
		 angle_3_object)
CGAL_Kernel_pred(CGALi::Collinear,
	         Collinear_3,
		 collinear_3_object)
CGAL_Kernel_pred(CGALi::Coplanar,
	         Coplanar_3,
		 coplanar_3_object)
CGAL_Kernel_pred(CGALi::Coplanar_orientation,
	         Coplanar_orientation_3,
		 coplanar_orientation_3_object)
CGAL_Kernel_pred(CGALi::Coplanar_side_of_bounded_circle,
	         Coplanar_side_of_bounded_circle_3,
		 coplanar_side_of_bounded_circle_3_object)
CGAL_Kernel_pred(CGAL ::p_Orientation<Point_3>,
	         Orientation_3,
		 orientation_3_object)
CGAL_Kernel_pred(CGALi::Call_is_degenerate,
	         Is_degenerate_3,
		 is_degenerate_3_object)
CGAL_Kernel_pred(CGALi::Call_has_on,
	         Has_on_3,
		 has_on_3_object)
CGAL_Kernel_pred(CGALi::Call_has_on_bounded_side,
	         Has_on_bounded_side_3,
		 has_on_bounded_side_3_object)
CGAL_Kernel_pred(CGALi::Call_has_on_unbounded_side,
	         Has_on_unbounded_side_3,
		 has_on_unbounded_side_3_object)
CGAL_Kernel_pred(CGALi::Call_has_on_boundary,
	         Has_on_boundary_3,
		 has_on_boundary_3_object)
CGAL_Kernel_pred(CGALi::Call_has_on_positive_side,
	         Has_on_positive_side_3,
		 has_on_positive_side_3_object)
CGAL_Kernel_pred(CGALi::Call_has_on_negative_side,
	         Has_on_negative_side_3,
		 has_on_negative_side_3_object)
CGAL_Kernel_pred(CGALi::Call_oriented_side,
	         Oriented_side_3,
		 oriented_side_3_object)
CGAL_Kernel_pred(CGALi::Call_bounded_side,
	         Bounded_side_3,
		 bounded_side_3_object)
CGAL_Kernel_pred(CGALi::Are_ordered_along_line,
	         Are_ordered_along_line_3,
		 are_ordered_along_line_3_object)
CGAL_Kernel_pred(CGALi::Are_strictly_ordered_along_line,
	         Are_strictly_ordered_along_line_3,
		 are_strictly_ordered_along_line_3_object)
CGAL_Kernel_pred(CGALi::Collinear_are_ordered_along_line,
	         Collinear_are_ordered_along_line_3,
		 collinear_are_ordered_along_line_3_object)
CGAL_Kernel_pred(CGALi::Collinear_are_strictly_ordered_along_line,
	         Collinear_are_strictly_ordered_along_line_3,
		 collinear_are_strictly_ordered_along_line_3_object)
CGAL_Kernel_pred(CGALi::Side_of_oriented_sphere,
	         Side_of_oriented_sphere_3,
		 side_of_oriented_sphere_3_object)
CGAL_Kernel_pred(CGALi::Side_of_bounded_sphere,
	         Side_of_bounded_sphere_3,
		 side_of_bounded_sphere_3_object)

CGAL_Kernel_cons(CGALi::Construct<Point_d>,
	         Construct_point_d,
		 construct_point_d_object)

#undef CGAL_Kernel_pred
#undef CGAL_Kernel_cons
#undef CGAL_Kernel_pred2
#undef CGAL_Kernel_cons2


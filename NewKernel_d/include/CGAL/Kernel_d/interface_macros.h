#ifndef CGAL_Kernel_pred
#  define CGAL_Kernel_pred(X, Y)
#endif

#ifndef CGAL_Kernel_comp
#  define CGAL_Kernel_comp(X, Y)
#endif
#ifndef CGAL_Kernel_comp1
#  define CGAL_Kernel_comp1(X, Y) CGAL_Kernel_comp(X, Y)
#endif
#ifndef CGAL_Kernel_comp2
#  define CGAL_Kernel_comp2(X, Y) CGAL_Kernel_comp(X, Y)
#endif

#ifndef CGAL_Kernel_cons
#  define CGAL_Kernel_cons(X, Y)
#endif
#ifndef CGAL_Kernel_cons1
#  define CGAL_Kernel_cons1(X, Y) CGAL_Kernel_cons(X, Y)
#endif
#ifndef CGAL_Kernel_cons2
#  define CGAL_Kernel_cons2(X, Y) CGAL_Kernel_cons(X, Y)
#endif

#ifndef CGAL_Kernel_obj
#  define CGAL_Kernel_obj(X,Y)
#endif
#ifndef CGAL_Kernel_obj1
#  define CGAL_Kernel_obj1(X,Y) CGAL_Kernel_obj(X,Y)
#endif
#ifndef CGAL_Kernel_obj2
#  define CGAL_Kernel_obj2(X,Y) CGAL_Kernel_obj(X,Y)
#endif

CGAL_Kernel_obj1(Vector,vector)
CGAL_Kernel_obj1(Point,point)
//CGAL_Kernel_obj2(Segment,segment)

CGAL_Kernel_cons1(Construct_point_cartesian_const_iterator,
		  construct_point_cartesian_const_iterator_object)
CGAL_Kernel_cons1(Construct_vector_cartesian_const_iterator,
		  construct_vector_cartesian_const_iterator_object)
CGAL_Kernel_cons2(Sum_of_vectors,
		  sum_of_vectors_object)
CGAL_Kernel_cons2(Difference_of_vectors,
		  difference_of_vectors_object)
CGAL_Kernel_cons2(Opposite_vector,
		  opposite_vector_object)
CGAL_Kernel_cons2(Difference_of_points,
		  difference_of_points_object)
CGAL_Kernel_cons2(Vector_to_point,
		  vector_to_point_object)
CGAL_Kernel_cons2(Point_to_vector,
		  point_to_vector_object)
CGAL_Kernel_cons2(Midpoint,
		  midpoint_object)
//CGAL_Kernel_cons2(Construct_segment,
//		  construct_segment_object)
//CGAL_Kernel_cons2(Segment_extremity,
//		  segment_extremity_object)

CGAL_Kernel_comp1(Compute_point_cartesian_coordinate,
		  compute_point_cartesian_coordinate_object)
CGAL_Kernel_comp1(Compute_vector_cartesian_coordinate,
		  compute_vector_cartesian_coordinate_object)
CGAL_Kernel_comp2(Squared_length,
		  squared_length_object)
CGAL_Kernel_comp2(Squared_distance,
		  squared_distance_object)
CGAL_Kernel_comp2(Squared_distance_to_origin,
		  squared_distance_to_origin_object)

#if 0
CGAL_Kernel_pred(Orientation,
		 orientation_object)
#endif
CGAL_Kernel_pred(Orientation_of_points,
		 orientation_of_points_object)
CGAL_Kernel_pred(Orientation_of_vectors,
		 orientation_of_vectors_object)
CGAL_Kernel_pred(Side_of_oriented_sphere,
		 side_of_oriented_sphere_object)
CGAL_Kernel_pred(Less_point_cartesian_coordinate,
		 less_point_cartesian_coordinate_object)
CGAL_Kernel_pred(Compare_point_cartesian_coordinate,
		 compare_point_cartesian_coordinate_object)
CGAL_Kernel_pred(Compare_distance,
		 compare_distance_object)
CGAL_Kernel_pred(Compare_lexicographically,
		 compare_lexicographically_object)
CGAL_Kernel_pred(Less_lexicographically,
		 less_lexicographically_object)
CGAL_Kernel_pred(Equal_points,
		 equal_points_object)
CGAL_Kernel_pred(Less_or_equal_lexicographically,
		 less_or_equal_lexicographically_object)
CGAL_Kernel_pred(Contained_in_affine_hull,
		 contained_in_affine_hull_object)
CGAL_Kernel_pred(In_flat_orientation,
		 in_flat_orientation_object)
CGAL_Kernel_pred(In_flat_side_of_oriented_sphere,
		 in_flat_side_of_oriented_sphere_object)
CGAL_Kernel_pred(Construct_flat_orientation,
		 construct_flat_orientation_object)

#undef CGAL_Kernel_pred
#undef CGAL_Kernel_comp
#undef CGAL_Kernel_comp1
#undef CGAL_Kernel_comp2
#undef CGAL_Kernel_cons
#undef CGAL_Kernel_cons1
#undef CGAL_Kernel_cons2
#undef CGAL_Kernel_obj
#undef CGAL_Kernel_obj1
#undef CGAL_Kernel_obj2

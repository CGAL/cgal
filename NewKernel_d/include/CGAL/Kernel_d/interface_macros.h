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
#  define CGAL_Kernel_obj(X)
#endif
#ifndef CGAL_Kernel_obj1
#  define CGAL_Kernel_obj1(X) CGAL_Kernel_obj(X)
#endif
#ifndef CGAL_Kernel_obj2
#  define CGAL_Kernel_obj2(X) CGAL_Kernel_obj(X)
#endif

CGAL_Kernel_obj1(Vector)
CGAL_Kernel_obj1(Point)
CGAL_Kernel_obj2(Segment)

CGAL_Kernel_cons1(Construct_vector,
		  construct_vector_object)
CGAL_Kernel_cons1(Construct_point,
		  construct_point_object)
CGAL_Kernel_cons1(Construct_cartesian_const_iterator,
		  construct_cartesian_const_iterator_object)
CGAL_Kernel_cons2(Construct_sum_of_vectors,
		  construct_sum_of_vectors_object)
CGAL_Kernel_cons2(Construct_difference_of_vectors,
		  construct_difference_of_vectors_object)
CGAL_Kernel_cons2(Construct_opposite_vector,
		  construct_opposite_vector_object)
CGAL_Kernel_cons2(Construct_midpoint,
		  construct_midpoint_object)
CGAL_Kernel_cons2(Construct_segment,
		  construct_segment_object)
CGAL_Kernel_cons2(Construct_segment_extremity,
		  construct_segment_extremity_object)

CGAL_Kernel_comp1(Compute_cartesian_coordinate,
		  compute_cartesian_coordinate_object)

CGAL_Kernel_pred(Orientation,
		 orientation_object)

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

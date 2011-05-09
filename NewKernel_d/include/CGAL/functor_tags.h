#ifndef CGAL_FUNCTOR_TAGS_H
#define CGAL_FUNCTOR_TAGS_H
namespace CGAL {
	class Null_type {~Null_type();}; // no such object should be created

	struct Vector_tag {};
	struct Point_tag {};
	struct Segment_tag {};
	struct Line_tag {};
	struct Direction_tag {};
	struct Ray_tag {};
	struct Bbox_tag {};

	struct Construct_vector_tag {};
	struct Construct_point_tag {};
	struct Construct_segment_tag {};
	struct Construct_line_tag {};
	struct Construct_direction_tag {};
	struct Construct_ray_tag {};
	struct Construct_cartesian_const_iterator_tag {};
	struct Construct_midpoint_tag {};
	struct Construct_sum_of_vectors_tag {};
	struct Construct_difference_of_vectors_tag {};
	struct Construct_opposite_vector_tag {};

	struct Compute_cartesian_coordinate_tag {};
	struct Compute_homogeneous_coordinate_tag {};
	struct Compute_squared_distance_tag {};
	struct Compute_squared_length_tag {};

	struct Predicate_less_cartesian_coordinate_tag {};
	//FIXME: choose a convention
	struct Predicate_orientation_tag {};
	struct Orientation_tag {};
	struct Predicate_in_sphere_tag {};
}
#endif // CGAL_FUNCTOR_TAGS_H

#ifndef CGAL_FUNCTOR_TAGS_H
#define CGAL_FUNCTOR_TAGS_H
namespace CGAL {
	class Null_type {~Null_type();}; // no such object should be created

	template<class,class>struct map_kernel_obj{typedef Null_type type;};
#define DECL_OBJ(X) struct X##_tag {}; \
	template<class K>struct map_kernel_obj<K,X##_tag>{typedef typename K::X type;}
	DECL_OBJ(Vector);
	DECL_OBJ(Point);
	DECL_OBJ(Segment);
	DECL_OBJ(Line);
	DECL_OBJ(Direction);
	DECL_OBJ(Ray);
	DECL_OBJ(Bbox);
#undef DECL_OBJ

	struct Construct_point_cartesian_const_iterator_tag {};
	struct Construct_vector_cartesian_const_iterator_tag {};

	template<class T>struct map_result_tag{typedef Null_type type;};
#define DECL_CONSTRUCT(X,Y) struct X##_tag {}; \
	template<>struct map_result_tag<X##_tag>{typedef Y##_tag type;}
	DECL_CONSTRUCT(Construct_vector,Vector);
	DECL_CONSTRUCT(Construct_point,Point);
	DECL_CONSTRUCT(Construct_segment,Segment);
	DECL_CONSTRUCT(Construct_line,Line);
	DECL_CONSTRUCT(Construct_direction,Direction);
	DECL_CONSTRUCT(Construct_ray,Ray);
	DECL_CONSTRUCT(Construct_midpoint,Point);
	DECL_CONSTRUCT(Construct_segment_extremity,Point);
	DECL_CONSTRUCT(Construct_sum_of_vectors,Vector);
	DECL_CONSTRUCT(Construct_difference_of_vectors,Vector);
	DECL_CONSTRUCT(Construct_opposite_vector,Vector);
#undef DECL_CONSTRUCT

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

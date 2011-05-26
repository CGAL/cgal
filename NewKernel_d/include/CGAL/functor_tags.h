#ifndef CGAL_FUNCTOR_TAGS_H
#define CGAL_FUNCTOR_TAGS_H
namespace CGAL {
	class Null_type {~Null_type();}; // no such object should be created

	struct Begin_tag {};
	struct End_tag {};

	struct Predicate_tag {};
	struct Construct_tag {};
	struct Compute_tag {};
	struct Misc_tag {};

	template<class>struct map_functor_type{typedef Misc_tag type;};

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

	template<class T>struct map_result_tag{typedef Null_type type;};
#define DECL_CONSTRUCT(X,Y) struct X##_tag {}; \
	template<>struct map_result_tag<X##_tag>{typedef Y##_tag type;}; \
	template<>struct map_functor_type<X##_tag>{typedef Construct_tag type;}
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

#define DECL_COMPUTE(X) struct X##_tag {}; \
	template<>struct map_functor_type<X##_tag>{typedef Compute_tag type;}
	DECL_COMPUTE(Compute_cartesian_coordinate);
	DECL_COMPUTE(Compute_homogeneous_coordinate);
	DECL_COMPUTE(Compute_squared_distance);
	DECL_COMPUTE(Compute_squared_length);
#undef DECL_COMPUTE

	//FIXME: choose a convention: prefix with Predicate_ ?
#define DECL_PREDICATE(X) struct X##_tag {}; \
	template<>struct map_functor_type<X##_tag>{typedef Predicate_tag type;}
	DECL_PREDICATE(Less_cartesian_coordinate);
	DECL_PREDICATE(Orientation);
	DECL_PREDICATE(Orientation_of_points);
	DECL_PREDICATE(Orientation_of_vectors);
	DECL_PREDICATE(In_sphere);
#undef DECL_PREDICATE

#define DECL_MISC(X) struct X##_tag {}; \
	template<>struct map_functor_type<X##_tag>{typedef Misc_tag type;}
	//TODO: split into _begin and _end ?
	DECL_MISC(Construct_point_cartesian_const_iterator);
	DECL_MISC(Construct_vector_cartesian_const_iterator);
	DECL_MISC(Point_dimension);
	DECL_MISC(Vector_dimension);
#undef DECL_MISC
}
#endif // CGAL_FUNCTOR_TAGS_H

#ifndef CGAL_FUNCTOR_TAGS_H
#define CGAL_FUNCTOR_TAGS_H
#include <boost/mpl/has_xxx.hpp>
namespace CGAL {
	class Null_type {~Null_type();}; // no such object should be created

	// To construct iterators
	struct Begin_tag {};
	struct End_tag {};

	// Functor category
	struct Predicate_tag {};
	struct Construct_tag {};
	struct Construct_iterator_tag {};
	struct Compute_tag {};
	struct Misc_tag {};

	struct No_filter_tag {};

	template<class>struct Construct_ttag {};
	template<class>struct Convert_ttag {};

	template<class> struct map_functor_type { typedef Misc_tag type; };
	template<class Tag, class Obj, class Base> struct Typedef_tag_type;
	template<class Kernel, class Tag> struct Read_tag_type {};
	template<class Kernel, class Tag> struct Provides_tag_type;


#define DECL_OBJ(X) struct X##_tag {}; \
  template<class Obj,class Base> \
  struct Typedef_tag_type<X##_tag, Obj, Base> : Base { typedef Obj X; }; \
  namespace has_object { BOOST_MPL_HAS_XXX_TRAIT_DEF(X) } \
  template<class Kernel> \
  struct Provides_tag_type<Kernel, X##_tag> : has_object::has_##X<Kernel> {}; \
  template<class Kernel> \
  struct Read_tag_type<Kernel, X##_tag> { typedef typename Kernel::X type; }

	// Not exactly objects, but the extras can't hurt.
	DECL_OBJ(FT);
	DECL_OBJ(RT);

	DECL_OBJ(Vector);
	DECL_OBJ(Point);
	DECL_OBJ(Segment);
	DECL_OBJ(Line);
	DECL_OBJ(Direction);
	DECL_OBJ(Ray);
	DECL_OBJ(Bbox);
#undef DECL_OBJ

	template<class> struct is_NT_tag { enum { value = false }; };
	template<> struct is_NT_tag<FT_tag> { enum { value = true }; };
	template<> struct is_NT_tag<RT_tag> { enum { value = true }; };

	template<class> struct iterator_tag_traits {
	  enum { is_iterator = false, has_nth_element = false };
	  typedef Null_tag value_tag;
	};

#define DECL_COMPUTE(X) struct X##_tag {}; \
	template<>struct map_functor_type<X##_tag>{typedef Compute_tag type;}
	DECL_COMPUTE(Compute_point_cartesian_coordinate);
	DECL_COMPUTE(Compute_vector_cartesian_coordinate);
	DECL_COMPUTE(Compute_homogeneous_coordinate);
	DECL_COMPUTE(Compute_squared_distance);
	DECL_COMPUTE(Compute_squared_length);
	DECL_COMPUTE(Compute_scalar_product);
#undef DECL_COMPUTE

#define DECL_ITER_OBJ(X,Y,Z,C) struct X##_tag {}; \
  template<>struct iterator_tag_traits<X##_tag> { \
    enum { is_iterator = true, has_nth_element = true }; \
    typedef Y##_tag value_tag; \
    typedef Z##_tag nth_element; \
    typedef C##_tag container; \
  }; \
  template<class Obj,class Base> \
  struct Typedef_tag_type<X##_tag, Obj, Base> : Base { typedef Obj X; }; \
  namespace has_object { BOOST_MPL_HAS_XXX_TRAIT_DEF(X) } \
  template<class Kernel> \
  struct Provides_tag_type<Kernel, X##_tag> : has_object::has_##X<Kernel> {}; \
  template<class Kernel> \
  struct Read_tag_type<Kernel, X##_tag> { typedef typename Kernel::X type; }

	DECL_ITER_OBJ(Vector_cartesian_const_iterator, FT, Compute_vector_cartesian_coordinate, Vector);
	DECL_ITER_OBJ(Point_cartesian_const_iterator, FT, Compute_point_cartesian_coordinate, Point);
#undef DECL_ITER_OBJ

	template<class>struct map_result_tag{typedef Null_type type;};
	template<class T>struct map_result_tag<Construct_ttag<T> >{typedef T type;};

	template<class T>struct map_functor_type<Construct_ttag<T> > :
	  BOOSTD conditional<iterator_tag_traits<T>::is_iterator,
		 Construct_iterator_tag,
		 Construct_tag> {};

	template<class T>struct map_functor_type<Convert_ttag<T> >{typedef Misc_tag type;};
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
#if 0
#define DECL_ITER_CONSTRUCT(X,Y) struct X##_tag {}; \
	template<>struct map_result_tag<X##_tag>{typedef Y##_tag type;}; \
	template<>struct map_functor_type<X##_tag>{typedef Construct_iterator_tag type;}
	DECL_ITER_CONSTRUCT(Construct_point_cartesian_const_iterator,Point_cartesian_const_iterator);
	DECL_ITER_CONSTRUCT(Construct_vector_cartesian_const_iterator,Vector_cartesian_const_iterator);
#undef DECL_ITER_CONSTRUCT
#endif

	//FIXME: choose a convention: prefix with Predicate_ ?
#define DECL_PREDICATE(X) struct X##_tag {}; \
	template<>struct map_functor_type<X##_tag>{typedef Predicate_tag type;}
	DECL_PREDICATE(Less_point_cartesian_coordinate);
	DECL_PREDICATE(Compare_point_cartesian_coordinate);
	DECL_PREDICATE(Compare_distance);
	DECL_PREDICATE(Compare_lexicographically);
	DECL_PREDICATE(Orientation);
	DECL_PREDICATE(Orientation_of_points);
	DECL_PREDICATE(Orientation_of_vectors);
	DECL_PREDICATE(Side_of_oriented_sphere);
	DECL_PREDICATE(Contained_in_affine_hull);
	DECL_PREDICATE(In_flat_orientation);
	DECL_PREDICATE(In_flat_side_of_oriented_sphere);
	DECL_PREDICATE(Construct_flat_orientation); // Making it a predicate is a questionable choice, it should be possible to let it be a construction for some implementations. Not sure how to do that... TODO
#undef DECL_PREDICATE

#define DECL_MISC(X) struct X##_tag {}; \
	template<>struct map_functor_type<X##_tag>{typedef Misc_tag type;}
	//TODO: split into _begin and _end ?
	//DECL_MISC(Construct_point_cartesian_const_iterator);
	//DECL_MISC(Construct_vector_cartesian_const_iterator);
	DECL_MISC(Point_dimension);
	DECL_MISC(Vector_dimension);
#undef DECL_MISC


	// Properties for LA
	struct Has_extra_dimension_tag {};
	struct Has_vector_plus_minus_tag {};
	struct Has_vector_scalar_ops_tag {};
	struct Has_dot_product_tag {};
	struct Has_determinant_of_vectors_tag {};
	struct Has_determinant_of_points_tag {};
	struct Has_determinant_of_iterator_to_vectors_tag {};
	struct Has_determinant_of_iterator_to_points_tag {};
	struct Has_determinant_of_vectors_omit_last_tag {};

	template<class> struct Preserved_by_non_linear_extra_coordinate
	  : boost::false_type {};
	template<> struct Preserved_by_non_linear_extra_coordinate
	  <Has_extra_dimension_tag> : boost::true_type {};
	template<> struct Preserved_by_non_linear_extra_coordinate
	  <Has_determinant_of_vectors_tag> : boost::true_type {};
	template<> struct Preserved_by_non_linear_extra_coordinate
	  <Has_determinant_of_points_tag> : boost::true_type {};
	template<> struct Preserved_by_non_linear_extra_coordinate
	  <Has_determinant_of_iterator_to_vectors_tag> : boost::true_type {};
	template<> struct Preserved_by_non_linear_extra_coordinate
	  <Has_determinant_of_iterator_to_points_tag> : boost::true_type {};
	template<> struct Preserved_by_non_linear_extra_coordinate
	  <Has_determinant_of_vectors_omit_last_tag> : boost::true_type {};

	// Kernel properties
	struct Point_stores_squared_distance_to_origin_tag {};

}
#endif // CGAL_FUNCTOR_TAGS_H

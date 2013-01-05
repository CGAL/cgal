#ifndef CGAL_FUNCTOR_TAGS_H
#define CGAL_FUNCTOR_TAGS_H
#include <CGAL/tags.h> // for Null_tag
#include <CGAL/marcutils.h>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/not.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>
namespace CGAL {

  // Find a better place for this later

  template <class K, class T, class=void> struct Get_type
    : K::template Type<T> {};
  template <class K, class F, class O=void, class=void> struct Get_functor
    : K::template Functor<F, O> {};
#ifdef CGAL_CXX0X
  template <class K, class T> using Type = typename Get_type<K, T>::type;
  template <class K, class T> using Functor = typename Get_functor<K, T>::type;
#endif

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

	template <class K, class F, class=void, class=void> struct Get_functor_category { typedef Misc_tag type; };
	template<class Tg, class Obj, class Base> struct Typedef_tag_type;
	//template<class Kernel, class Tg> struct Read_tag_type {};

	template<class Kernel, class Tg>
	struct Provides_type
	  : Has_type_different_from<Get_type<Kernel, Tg>, Null_type> {};

	template<class Kernel, class Tg, class O=void>
	struct Provides_functor
	  : Has_type_different_from<Get_functor<Kernel, Tg, O>, Null_functor> {};

	template<class K, class List, bool=boost::mpl::empty<List>::type::value>
	struct Provides_functors : boost::mpl::and_ <
				   Provides_functor<K, typename boost::mpl::front<List>::type>,
				   Provides_functors<K, typename boost::mpl::pop_front<List>::type> > {};
	template<class K, class List>
	struct Provides_functors<K, List, true> : boost::true_type {};

	template<class K, class List, bool=boost::mpl::empty<List>::type::value>
	struct Provides_types : boost::mpl::and_ <
				   Provides_type<K, typename boost::mpl::front<List>::type>,
				   Provides_types<K, typename boost::mpl::pop_front<List>::type> > {};
	template<class K, class List>
	struct Provides_types<K, List, true> : boost::true_type {};

	namespace internal { BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_Type,template Type<Null_tag>,false) }
	template<class Kernel, class Tg,
	  bool = internal::has_Type<Kernel>::value /* false */>
	struct Provides_type_i : boost::false_type {};
	template<class Kernel, class Tg>
	struct Provides_type_i <Kernel, Tg, true>
	  : Has_type_different_from<typename Kernel::template Type<Tg>, Null_type> {};

	namespace internal { BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_Functor,template Functor<Null_tag>,false) }
	template<class Kernel, class Tg, class O=void,
	  bool = internal::has_Functor<Kernel>::value /* false */>
	struct Provides_functor_i : boost::false_type {};
	template<class Kernel, class Tg, class O>
	struct Provides_functor_i <Kernel, Tg, O, true>
	  : Has_type_different_from<typename Kernel::template Functor<Tg, O>, Null_functor> {};

	struct Number_tag {};
	struct Discrete_tag {};
	struct Object_tag {};
	template <class K, class T, class=void> struct Get_type_category {
	  // The lazy kernel uses it too eagerly,
	  // so it currently needs a default.
	  typedef Null_tag type;
	};

#define DECL_OBJ_(X,Y) \
  template<class Obj,class Base> \
  struct Typedef_tag_type<X##_tag, Obj, Base> : Base { typedef Obj X; }; \
  template<class K, class D> \
  struct Get_type_category <K, X##_tag, D> { typedef Y##_tag type; }
#define DECL_OBJ(X,Y) struct X##_tag {}; \
  DECL_OBJ_(X,Y)

  //namespace has_object { BOOST_MPL_HAS_XXX_TRAIT_DEF(X) }
  //template<class Kernel>
  //struct Provides_tag_type<Kernel, X##_tag> : has_object::has_##X<Kernel> {};
  //template<class Kernel>
  //struct Read_tag_type<Kernel, X##_tag> { typedef typename Kernel::X type; }

	// Not exactly objects, but the extras can't hurt.
	DECL_OBJ(FT, Number);
	DECL_OBJ(RT, Number);

	DECL_OBJ(Bool, Discrete); // Boolean_tag is already taken, and is a template :-(
	DECL_OBJ(Comparison_result, Discrete);
	DECL_OBJ(Sign, Discrete);
	DECL_OBJ(Orientation, Discrete); // Note: duplicate with the functor tag!
	DECL_OBJ(Oriented_side, Discrete);
	DECL_OBJ(Bounded_side, Discrete);
	DECL_OBJ(Angle, Discrete);
	DECL_OBJ(Flat_orientation, Discrete);

	DECL_OBJ(Vector, Object);
	DECL_OBJ(Point, Object);
	DECL_OBJ(Segment, Object);
	DECL_OBJ(Sphere, Object);
	DECL_OBJ(Line, Object);
	DECL_OBJ(Direction, Object);
	DECL_OBJ(Hyperplane, Object);
	DECL_OBJ(Ray, Object);
	DECL_OBJ(Iso_box, Object);
	DECL_OBJ(Bbox, Object);
	DECL_OBJ(Aff_transformation, Object);
#undef DECL_OBJ_
#undef DECL_OBJ

	CGAL_KD_DEFAULT_TYPE(RT_tag,(typename Get_type<K, FT_tag>::type),(),());
	CGAL_KD_DEFAULT_TYPE(FT_tag,(CGAL::Quotient<typename Get_type<K, RT_tag>::type>),(),());

#define SMURF2(A,B) CGAL_KD_DEFAULT_TYPE(A##_tag,(typename Same_uncertainty_nt<B, typename Get_type<K,RT_tag>::type>::type),(RT_tag),())
#define SMURF1(A) SMURF2(A,CGAL::A)
	SMURF2(Bool, bool);
	SMURF1(Sign);
	SMURF1(Comparison_result);
	SMURF1(Orientation);
	SMURF1(Oriented_side);
	SMURF1(Bounded_side);
	SMURF1(Angle);
#undef SMURF1
#undef SMURF2

	// TODO: replace with Get_type_category
	template<class> struct is_NT_tag { enum { value = false }; };
	template<> struct is_NT_tag<FT_tag> { enum { value = true }; };
	template<> struct is_NT_tag<RT_tag> { enum { value = true }; };

	template<class> struct iterator_tag_traits {
	  enum { is_iterator = false, has_nth_element = false };
	  typedef Null_tag value_tag;
	};

#define DECL_COMPUTE(X) struct X##_tag {}; \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Compute_tag type;}
	DECL_COMPUTE(Compute_point_cartesian_coordinate);
	DECL_COMPUTE(Compute_vector_cartesian_coordinate);
	DECL_COMPUTE(Compute_homogeneous_coordinate);
	DECL_COMPUTE(Squared_distance);
	DECL_COMPUTE(Squared_distance_to_origin);
	DECL_COMPUTE(Squared_length);
	DECL_COMPUTE(Squared_radius);
	DECL_COMPUTE(Scalar_product);
	DECL_COMPUTE(Hyperplane_translation);
#undef DECL_COMPUTE

#define DECL_ITER_OBJ(X,Y,Z,C) struct X##_tag {}; \
  template<>struct iterator_tag_traits<X##_tag> { \
    enum { is_iterator = true, has_nth_element = true }; \
    typedef Y##_tag value_tag; \
    typedef Z##_tag nth_element; \
    typedef C##_tag container; \
  }; \
  template<class Obj,class Base> \
  struct Typedef_tag_type<X##_tag, Obj, Base> : Base { typedef Obj X; }

  //namespace has_object { BOOST_MPL_HAS_XXX_TRAIT_DEF(X) }
  //template<class Kernel>
  //struct Provides_tag_type<Kernel, X##_tag> : has_object::has_##X<Kernel> {};
  //template<class Kernel>
  //struct Read_tag_type<Kernel, X##_tag> { typedef typename Kernel::X type; }

	DECL_ITER_OBJ(Vector_cartesian_const_iterator, FT, Compute_vector_cartesian_coordinate, Vector);
	DECL_ITER_OBJ(Point_cartesian_const_iterator, FT, Compute_point_cartesian_coordinate, Point);
#undef DECL_ITER_OBJ

	template<class>struct map_result_tag{typedef Null_type type;};
	template<class T>struct map_result_tag<Construct_ttag<T> >{typedef T type;};

	template<class A,class T,class B,class C>struct Get_functor_category<A,Construct_ttag<T>,B,C> :
	  BOOSTD conditional<iterator_tag_traits<T>::is_iterator,
		 Construct_iterator_tag,
		 Construct_tag> {};

	// Really?
	template<class A,class T,class B,class C>struct Get_functor_category<A,Convert_ttag<T>,B,C>{typedef Misc_tag type;};
#define DECL_CONSTRUCT(X,Y) struct X##_tag {}; \
	template<>struct map_result_tag<X##_tag>{typedef Y##_tag type;}; \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Construct_tag type;}
	DECL_CONSTRUCT(Midpoint,Point);
	DECL_CONSTRUCT(Center_of_sphere,Point);
	DECL_CONSTRUCT(Segment_extremity,Point);
	DECL_CONSTRUCT(Sum_of_vectors,Vector);
	DECL_CONSTRUCT(Difference_of_vectors,Vector);
	DECL_CONSTRUCT(Opposite_vector,Vector);
	DECL_CONSTRUCT(Scaled_vector,Vector);
	DECL_CONSTRUCT(Orthogonal_vector,Vector);
	DECL_CONSTRUCT(Difference_of_points,Vector);
	DECL_CONSTRUCT(Translated_point,Point);
	DECL_CONSTRUCT(Point_to_vector,Vector);
	DECL_CONSTRUCT(Vector_to_point,Point);
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
#define DECL_PREDICATE_(X) \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Predicate_tag type;}
#define DECL_PREDICATE(X) struct X##_tag {}; \
	DECL_PREDICATE_(X)
	DECL_PREDICATE(Less_point_cartesian_coordinate);
	DECL_PREDICATE(Compare_point_cartesian_coordinate);
	DECL_PREDICATE(Compare_distance);
	DECL_PREDICATE(Compare_lexicographically);
	DECL_PREDICATE(Less_lexicographically);
	DECL_PREDICATE(Less_or_equal_lexicographically);
	DECL_PREDICATE(Equal_points);
	DECL_PREDICATE_(Orientation); // duplicate with the type
	DECL_PREDICATE(Orientation_of_points);
	DECL_PREDICATE(Orientation_of_vectors);
	DECL_PREDICATE(Side_of_oriented_sphere);
	DECL_PREDICATE(Contained_in_affine_hull);
	DECL_PREDICATE(In_flat_orientation);
	DECL_PREDICATE(In_flat_side_of_oriented_sphere);
	DECL_PREDICATE(Construct_flat_orientation); // Making it a predicate is a questionable choice, it should be possible to let it be a construction for some implementations. Not sure how to do that... TODO
#undef DECL_PREDICATE

#define DECL_MISC(X) struct X##_tag {}; \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Misc_tag type;}
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
	struct Stores_squared_norm_tag {};

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

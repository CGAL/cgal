// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Marc Glisse

#ifndef CGAL_FUNCTOR_TAGS_H
#define CGAL_FUNCTOR_TAGS_H
#include <CGAL/tags.h> // for Null_tag
#include <CGAL/NewKernel_d/utils.h>
#ifdef CGAL_CXX11
#include <type_traits>
#include <utility>
#endif
#include <boost/type_traits.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/if.hpp>
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
#ifdef CGAL_CXX11
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

	//// This version does not like Functor<T,bool=false>
	//namespace internal { BOOST_MPL_HAS_XXX_TEMPLATE_NAMED_DEF(has_Functor,Functor,false) }
	// This version lets us use non-type template parameters, but fails with older EDG-based compilers (Intel 14).
	namespace internal { BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_Functor,template Functor<Null_tag>,false) }

	template<class Kernel, class Tg, class O=void,
	  bool = internal::has_Functor<Kernel>::value /* false */>
	struct Provides_functor_i : boost::false_type {};
	template<class Kernel, class Tg, class O>
	struct Provides_functor_i <Kernel, Tg, O, true>
	  : Has_type_different_from<typename Kernel::template Functor<Tg, O>, Null_functor> {};

	// TODO: Refine this a bit.
	template <class K, class T, class D=void,
		  //bool=Provides_functor<K,T>::value,
		  //bool=Provides_functor_i<K,T>::value,
		  bool = internal::has_Functor<K>::value
		  >
	struct Inherit_functor : K::template Functor<T> {};
	template <class K, class T, class D>
	struct Inherit_functor <K, T, D, false> {};

	template <class K, class T, bool=internal::has_Type<K>::value>
	struct Inherit_type : K::template Type<T> {};
	template <class K, class T>
	struct Inherit_type <K, T, false> {};

	struct Number_tag {};
	struct Discrete_tag {};
	struct Object_tag {};
	template <class K, class T, class=void> struct Get_type_category {
	  // The lazy kernel uses it too eagerly,
	  // so it currently needs a default.
	  typedef Null_tag type;
	};

#define CGAL_DECL_OBJ_(X,Y) \
  template<class Obj,class Base> \
  struct Typedef_tag_type<X##_tag, Obj, Base> : Base { typedef Obj X; }; \
  template<class K, class D> \
  struct Get_type_category <K, X##_tag, D> { typedef Y##_tag type; }
#define CGAL_DECL_OBJ(X,Y) struct X##_tag {}; \
  CGAL_DECL_OBJ_(X,Y)

  //namespace has_object { BOOST_MPL_HAS_XXX_TRAIT_DEF(X) }
  //template<class Kernel>
  //struct Provides_tag_type<Kernel, X##_tag> : has_object::has_##X<Kernel> {};
  //template<class Kernel>
  //struct Read_tag_type<Kernel, X##_tag> { typedef typename Kernel::X type; }

	// Not exactly objects, but the extras can't hurt.
	CGAL_DECL_OBJ(FT, Number);
	CGAL_DECL_OBJ(RT, Number);

	CGAL_DECL_OBJ(Bool, Discrete); // Boolean_tag is already taken, and is a template :-(
	CGAL_DECL_OBJ(Comparison_result, Discrete);
	CGAL_DECL_OBJ(Sign, Discrete);
	CGAL_DECL_OBJ(Orientation, Discrete); // Note: duplicate with the functor tag!
	CGAL_DECL_OBJ(Oriented_side, Discrete);
	CGAL_DECL_OBJ(Bounded_side, Discrete);
	CGAL_DECL_OBJ(Angle, Discrete);
	CGAL_DECL_OBJ(Flat_orientation, Discrete);

	CGAL_DECL_OBJ(Vector, Object);
	CGAL_DECL_OBJ(Point, Object);
	CGAL_DECL_OBJ(Segment, Object);
	CGAL_DECL_OBJ(Sphere, Object);
	CGAL_DECL_OBJ(Line, Object);
	CGAL_DECL_OBJ(Direction, Object);
	CGAL_DECL_OBJ(Hyperplane, Object);
	CGAL_DECL_OBJ(Ray, Object);
	CGAL_DECL_OBJ(Iso_box, Object);
	CGAL_DECL_OBJ(Bbox, Object);
	CGAL_DECL_OBJ(Aff_transformation, Object);
	CGAL_DECL_OBJ(Weighted_point, Object);
#undef CGAL_DECL_OBJ_
#undef CGAL_DECL_OBJ

// Intel fails with those, and they are not so useful.
//	CGAL_KD_DEFAULT_TYPE(RT_tag,(typename Get_type<K, FT_tag>::type),(),());
//	CGAL_KD_DEFAULT_TYPE(FT_tag,(CGAL::Quotient<typename Get_type<K, RT_tag>::type>),(),());

#define CGAL_SMURF2(A,B) CGAL_KD_DEFAULT_TYPE(A##_tag,(typename Same_uncertainty_nt<B, typename Get_type<K,RT_tag>::type>::type),(RT_tag),())
#define CGAL_SMURF1(A) CGAL_SMURF2(A,CGAL::A)
	CGAL_SMURF2(Bool, bool);
	CGAL_SMURF1(Sign);
	CGAL_SMURF1(Comparison_result);
	CGAL_SMURF1(Orientation);
	CGAL_SMURF1(Oriented_side);
	CGAL_SMURF1(Bounded_side);
	CGAL_SMURF1(Angle);
#undef CGAL_SMURF1
#undef CGAL_SMURF2

	// TODO: replace with Get_type_category
	template<class> struct is_NT_tag { enum { value = false }; };
	template<> struct is_NT_tag<FT_tag> { enum { value = true }; };
	template<> struct is_NT_tag<RT_tag> { enum { value = true }; };

	template<class> struct iterator_tag_traits {
	  enum { is_iterator = false, has_nth_element = false };
	  typedef Null_tag value_tag;
	};

#define CGAL_DECL_COMPUTE(X) struct X##_tag {}; \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Compute_tag type;}
	CGAL_DECL_COMPUTE(Compute_point_cartesian_coordinate);
	CGAL_DECL_COMPUTE(Compute_vector_cartesian_coordinate);
	CGAL_DECL_COMPUTE(Compute_homogeneous_coordinate);
	CGAL_DECL_COMPUTE(Squared_distance);
	CGAL_DECL_COMPUTE(Squared_distance_to_origin);
	CGAL_DECL_COMPUTE(Squared_length);
	CGAL_DECL_COMPUTE(Squared_radius);
	CGAL_DECL_COMPUTE(Squared_circumradius);
	CGAL_DECL_COMPUTE(Scalar_product);
	CGAL_DECL_COMPUTE(Hyperplane_translation);
	CGAL_DECL_COMPUTE(Value_at);
	CGAL_DECL_COMPUTE(Point_weight);
	CGAL_DECL_COMPUTE(Power_distance);
	CGAL_DECL_COMPUTE(Power_distance_to_point);
#undef CGAL_DECL_COMPUTE

#define CGAL_DECL_ITER_OBJ(X,Y,Z,C) struct X##_tag {}; \
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

	CGAL_DECL_ITER_OBJ(Vector_cartesian_const_iterator, FT, Compute_vector_cartesian_coordinate, Vector);
	CGAL_DECL_ITER_OBJ(Point_cartesian_const_iterator, FT, Compute_point_cartesian_coordinate, Point);
#undef CGAL_DECL_ITER_OBJ

	template<class>struct map_result_tag{typedef Null_type type;};
	template<class T>struct map_result_tag<Construct_ttag<T> >{typedef T type;};

	template<class A,class T,class B,class C>struct Get_functor_category<A,Construct_ttag<T>,B,C> :
	  boost::mpl::if_c<iterator_tag_traits<T>::is_iterator,
			   Construct_iterator_tag,
			   Construct_tag> {};

	// Really?
	template<class A,class T,class B,class C>struct Get_functor_category<A,Convert_ttag<T>,B,C>{typedef Misc_tag type;};

#define CGAL_DECL_CONSTRUCT(X,Y) struct X##_tag {}; \
	template<>struct map_result_tag<X##_tag>{typedef Y##_tag type;}; \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Construct_tag type;}
	CGAL_DECL_CONSTRUCT(Midpoint,Point);
	CGAL_DECL_CONSTRUCT(Center_of_sphere,Point);
	CGAL_DECL_CONSTRUCT(Point_of_sphere,Point);
	CGAL_DECL_CONSTRUCT(Segment_extremity,Point);
	CGAL_DECL_CONSTRUCT(Sum_of_vectors,Vector);
	CGAL_DECL_CONSTRUCT(Difference_of_vectors,Vector);
	CGAL_DECL_CONSTRUCT(Opposite_vector,Vector);
	CGAL_DECL_CONSTRUCT(Scaled_vector,Vector);
	CGAL_DECL_CONSTRUCT(Orthogonal_vector,Vector);
	CGAL_DECL_CONSTRUCT(Difference_of_points,Vector);
	CGAL_DECL_CONSTRUCT(Translated_point,Point);
	CGAL_DECL_CONSTRUCT(Point_to_vector,Vector);
	CGAL_DECL_CONSTRUCT(Vector_to_point,Point);
	CGAL_DECL_CONSTRUCT(Construct_min_vertex,Point);
	CGAL_DECL_CONSTRUCT(Construct_max_vertex,Point);
	CGAL_DECL_CONSTRUCT(Construct_circumcenter,Point);
	CGAL_DECL_CONSTRUCT(Point_drop_weight,Point);
	CGAL_DECL_CONSTRUCT(Power_center,Weighted_point);
#undef CGAL_DECL_CONSTRUCT
#if 0
#define CGAL_DECL_ITER_CONSTRUCT(X,Y) struct X##_tag {}; \
	template<>struct map_result_tag<X##_tag>{typedef Y##_tag type;}; \
	template<>struct map_functor_type<X##_tag>{typedef Construct_iterator_tag type;}
	CGAL_DECL_ITER_CONSTRUCT(Construct_point_cartesian_const_iterator,Point_cartesian_const_iterator);
	CGAL_DECL_ITER_CONSTRUCT(Construct_vector_cartesian_const_iterator,Vector_cartesian_const_iterator);
#undef CGAL_DECL_ITER_CONSTRUCT
#endif

	//FIXME: choose a convention: prefix with Predicate_ ?
#define CGAL_DECL_PREDICATE_(X) \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Predicate_tag type;}
#define CGAL_DECL_PREDICATE(X) struct X##_tag {}; \
	CGAL_DECL_PREDICATE_(X)
	CGAL_DECL_PREDICATE(Less_point_cartesian_coordinate);
	CGAL_DECL_PREDICATE(Compare_point_cartesian_coordinate);
	CGAL_DECL_PREDICATE(Compare_distance);
	CGAL_DECL_PREDICATE(Compare_lexicographically);
	CGAL_DECL_PREDICATE(Less_lexicographically);
	CGAL_DECL_PREDICATE(Less_or_equal_lexicographically);
	CGAL_DECL_PREDICATE(Equal_points);
	CGAL_DECL_PREDICATE(Has_on_positive_side);
	CGAL_DECL_PREDICATE_(Orientation); // duplicate with the type
	CGAL_DECL_PREDICATE_(Oriented_side); // duplicate with the type
	CGAL_DECL_PREDICATE(Orientation_of_points);
	CGAL_DECL_PREDICATE(Orientation_of_vectors);
	CGAL_DECL_PREDICATE(Side_of_oriented_sphere);
	CGAL_DECL_PREDICATE(Side_of_bounded_sphere);
	CGAL_DECL_PREDICATE(Side_of_bounded_circumsphere);
	CGAL_DECL_PREDICATE(Contained_in_affine_hull);
	CGAL_DECL_PREDICATE(In_flat_orientation);
	CGAL_DECL_PREDICATE(In_flat_side_of_oriented_sphere);
	CGAL_DECL_PREDICATE(Construct_flat_orientation); // Making it a predicate is a questionable choice, it should be possible to let it be a construction for some implementations. Not sure how to do that... TODO
	CGAL_DECL_PREDICATE(Linear_rank);
	CGAL_DECL_PREDICATE(Affine_rank);
	CGAL_DECL_PREDICATE(Linearly_independent);
	CGAL_DECL_PREDICATE(Affinely_independent);
	CGAL_DECL_PREDICATE(Contained_in_linear_hull);
	CGAL_DECL_PREDICATE(Contained_in_simplex);
	CGAL_DECL_PREDICATE(Power_side_of_power_sphere_raw);
	CGAL_DECL_PREDICATE(Power_side_of_power_sphere);
	CGAL_DECL_PREDICATE(In_flat_power_side_of_power_sphere_raw);
	CGAL_DECL_PREDICATE(In_flat_power_side_of_power_sphere);
#undef CGAL_DECL_PREDICATE

#define CGAL_DECL_MISC(X) struct X##_tag {}; \
	template<class A,class B,class C>struct Get_functor_category<A,X##_tag,B,C>{typedef Misc_tag type;}
	//TODO: split into _begin and _end ?
	//CGAL_DECL_MISC(Construct_point_cartesian_const_iterator);
	//CGAL_DECL_MISC(Construct_vector_cartesian_const_iterator);
	CGAL_DECL_MISC(Point_dimension);
	CGAL_DECL_MISC(Vector_dimension);
	CGAL_DECL_MISC(Linear_base); // Find a more appropriate category?
#undef CGAL_DECL_MISC


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

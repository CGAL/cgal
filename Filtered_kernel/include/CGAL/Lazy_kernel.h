// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
//
// Author(s)     : Andreas Fabri, Sylvain Pion

#ifndef CGAL_LAZY_KERNEL_H
#define CGAL_LAZY_KERNEL_H

#include <CGAL/basic.h>
//#include <CGAL/Filtered_predicate.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_2.h>
#include <CGAL/Filtered_kernel/Cartesian_coordinate_iterator_3.h>
#include <CGAL/Lazy.h>

#include <boost/none.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4348) // redefinition of default parameter in nested template class
#endif

namespace CGAL {

namespace internal {

// SFINAE way to detect result_type typedefs.
template<typename T>
class Has_result_type_helper
{
  typedef char one;
  typedef struct { char arr[2]; } two;

  template<typename _Up>
  struct Wrapper {};

  template<typename U>
  static one test(Wrapper<typename U::result_type>*);

  template<typename U>
  static two test(...);

public:
  static const bool value = sizeof(test<T>(0)) == 1;
};

template<typename T>
struct Has_result_type 
  : boost::integral_constant< bool, 
                              Has_result_type_helper< typename boost::remove_cv<T>::type>::value>
{};

template <typename T>
struct Get_result_type {
  typedef typename T::result_type type;
};

template <typename T>
struct Lazy_result_type 
  : boost::mpl::eval_if< Has_result_type<T>, 
                         Get_result_type<T>, 
                         boost::mpl::identity<void> > 
{};

class Enum_holder {
protected:
  enum { NONE, NT, VARIANT, OBJECT, BBOX };
};

} // internal

// Exact_kernel = exact kernel that will be made lazy
// Kernel = lazy kernel

// the Generic base simplies applies the generic magic functor stupidly.
// then the real base fixes up a few special cases.
template < typename EK_, typename AK_, typename E2A_, typename Kernel_ >
class Lazy_kernel_generic_base : protected internal::Enum_holder
  // : public Filtered_kernel_base<EK_>
    // TODO : Static_filters_base too ?  Check performance
{
public:

  typedef AK_   Approximate_kernel;
  typedef EK_   Exact_kernel;
  typedef E2A_  E2A;
  typedef Kernel_ Kernel;

  // synonym identical to Filtered_kernel
  typedef AK_   FK;

  // Note: Approx_converter and Exact_converter are defined in <CGAL/Lazy.h>
  typedef Approx_converter<Kernel, Approximate_kernel>   C2F;
  typedef Exact_converter<Kernel, Exact_kernel>    C2E;

  template < typename Kernel2 >
  struct Base { typedef Lazy_kernel_generic_base<Exact_kernel, Approximate_kernel, E2A, Kernel2>  Type; };

  template < typename T >
  struct Ambient_dimension {
    typedef typename T::Ambient_dimension type;
  };

  template < typename T >
  struct Feature_dimension {
    typedef typename T::Feature_dimension type;
  };

  typedef typename Exact_kernel::Kernel_tag                       Kernel_tag;
  typedef typename Exact_kernel::Rep_tag                          Rep_tag;

  enum { Has_filtered_predicates = true };
  enum { Has_static_filters = false };
  typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

  // Types
  typedef CGAL::Lazy_exact_nt<typename Exact_kernel::FT>  FT;
  typedef FT RT;

  typedef typename Same_uncertainty_nt<bool, FT>::type
	                                                              Boolean;
  typedef typename Same_uncertainty_nt<CGAL::Sign, FT>::type
	                                                              Sign;
  typedef typename Same_uncertainty_nt<CGAL::Comparison_result, FT>::type
	                                                              Comparison_result;
  typedef typename Same_uncertainty_nt<CGAL::Orientation, FT>::type
		                                                      Orientation;
  typedef typename Same_uncertainty_nt<CGAL::Oriented_side, FT>::type
	                                                              Oriented_side;
  typedef typename Same_uncertainty_nt<CGAL::Bounded_side, FT>::type
	                                                              Bounded_side;
  typedef typename Same_uncertainty_nt<CGAL::Angle, FT>::type
	                                                              Angle;

  typedef CGAL::Object Object_2;
  typedef CGAL::Object Object_3;

#define CGAL_Kernel_obj(X) \
  typedef Lazy<typename Approximate_kernel::X, typename Exact_kernel::X, typename Exact_kernel::FT, E2A>  X;

  CGAL_Kernel_obj(Data_accessor_2)
  CGAL_Kernel_obj(Conic_2)

  typedef Cartesian_coordinate_iterator_2<Kernel> Cartesian_const_iterator_2;
  typedef Cartesian_coordinate_iterator_3<Kernel> Cartesian_const_iterator_3;

  // Aff_transformation_2/3 operations are not functorized, so treat it as
  // an exterior object for now.
  // CGAL_Kernel_obj(Aff_transformation_2)
  // CGAL_Kernel_obj(Aff_transformation_3)
  typedef CGAL::Aff_transformationC2<Kernel>              Aff_transformation_2;
  typedef CGAL::Aff_transformationC3<Kernel>              Aff_transformation_3;

private:
  // We use a combination of partial and logic to extract the right
  // construction. Constructions without a result_type always have to
  // be done through specializations.
  // 
  // The case distinction goes as follows: 
  // result_type == FT                              => NT
  // result_type == Object                          => Object
  // result_type == Bbox_2 || result_type == Bbox_3 => BBOX
  // default                                        => NONE
  // no result_type                                 => NONE
  // 
  //
  // we require a Dummy because we cannot have complete
  // specializations inside a non-namespace scope.
  // The default implementation does some default handling,
  // the special cases are filtered by partial specializations.
  template <typename Construction, typename Dummy = boost::none_t>
  struct Lazy_wrapper_traits : 
    boost::mpl::eval_if< internal::Has_result_type<Construction>,
                         boost::mpl::eval_if< boost::is_same< typename boost::remove_cv< 
                                                                typename boost::remove_reference< 
                                                                  typename internal::Lazy_result_type<Construction>::type 
                                                                  >::type >::type,
                                                              typename Approximate_kernel::FT>,
                                              boost::mpl::int_<NT>,
                                              boost::mpl::eval_if< boost::is_same< typename internal::Lazy_result_type<Construction>::type,
                                                                                   CGAL::Object >,
                                                                   boost::mpl::int_<OBJECT>,
                                                                   boost::mpl::eval_if< boost::mpl::or_< 
                                                                                          boost::is_same< typename internal::Lazy_result_type<Construction>::type, CGAL::Bbox_2 >, 
                                                                                          boost::is_same< typename internal::Lazy_result_type<Construction>::type, CGAL::Bbox_3 > >,
                                                                                        boost::mpl::int_<BBOX>,
                                                                                        boost::mpl::int_<NONE> > > >,
                         boost::mpl::int_<NONE> >::type {};
  
#define CGAL_WRAPPER_TRAIT(NAME, WRAPPER)                               \
  template<typename Dummy>                                              \
  struct Lazy_wrapper_traits<typename Approximate_kernel::NAME, Dummy>  \
    : boost::mpl::int_<WRAPPER> {};

#if CGAL_INTERSECTION_VERSION > 1
  CGAL_WRAPPER_TRAIT(Intersect_2, VARIANT)
  CGAL_WRAPPER_TRAIT(Intersect_3, VARIANT)
#else
  CGAL_WRAPPER_TRAIT(Intersect_2, OBJECT)
  CGAL_WRAPPER_TRAIT(Intersect_3, OBJECT)
#endif
  CGAL_WRAPPER_TRAIT(Compute_squared_radius_2, NT)
  CGAL_WRAPPER_TRAIT(Compute_x_3, NT)
  CGAL_WRAPPER_TRAIT(Compute_y_3, NT)
  CGAL_WRAPPER_TRAIT(Compute_z_3, NT)

#undef CGAL_WRAPPER_TRAIT

  template <typename Construction, int Type = Lazy_wrapper_traits<Construction>::value>
  struct Select_wrapper_impl;

  template <typename Construction>
  struct Select_wrapper_impl<Construction, NONE> {
    template<typename Kernel, typename AKC, typename EKC>
    struct apply { typedef Lazy_construction<Kernel, AKC, EKC> type; };
  };

  template <typename Construction>
  struct Select_wrapper_impl<Construction, NT> {
    template<typename Kernel, typename AKC, typename EKC>
    struct apply { typedef Lazy_construction_nt<Kernel, AKC, EKC> type; };
  };

  template <typename Construction>
  struct Select_wrapper_impl<Construction, VARIANT> {
    template<typename Kernel, typename AKC, typename EKC>
    struct apply { typedef Lazy_construction_variant<Kernel, AKC, EKC> type; };
  };

  template <typename Construction>
  struct Select_wrapper_impl<Construction, OBJECT> {
    template<typename Kernel, typename AKC, typename EKC>
    struct apply { typedef Lazy_construction_object<Kernel, AKC, EKC> type; };
  };

  template <typename Construction>
  struct Select_wrapper_impl<Construction, BBOX> {
    template<typename Kernel, typename AKC, typename EKC>
    struct apply { typedef Lazy_construction_bbox<Kernel, AKC, EKC> type; };
  };
  
  template <typename Construction>
  struct Select_wrapper : Select_wrapper_impl<Construction> {};

public:

  // We don't touch the predicates.
  // FIXME TODO : better use a layer of Filtered_kernel on top of everything,
  //              so that semi-static filters are used as well (?).
#define CGAL_Kernel_pred(P, Pf)  \
    typedef Filtered_predicate<typename Exact_kernel::P, typename Approximate_kernel::P, C2E, C2F> P; \
    P Pf() const { return P(); }

#define CGAL_Kernel_cons(C, Cf) \
  typedef typename Select_wrapper<typename Approximate_kernel::C>::template apply<Kernel, typename Approximate_kernel::C, typename Exact_kernel::C>::type C; \
  C Cf() const { return C(); }

#include <CGAL/Kernel/interface_macros.h>
};

template < typename EK_, typename AK_, typename E2A_, typename Kernel_ >
class Lazy_kernel_base
  : public Lazy_kernel_generic_base<EK_, AK_, E2A_, Kernel_>
{
public:
  typedef Kernel_ Kernel;
  typedef AK_   Approximate_kernel;
  typedef EK_   Exact_kernel;
  typedef E2A_  E2A;

  template < typename Kernel2 >
  struct Base { typedef Lazy_kernel_base<Exact_kernel, Approximate_kernel, E2A, Kernel2>  Type; };

  typedef CommonKernelFunctors::Assign_2<Kernel>        Assign_2;
  typedef CommonKernelFunctors::Assign_3<Kernel>        Assign_3;
  typedef Lazy_construction_bbox<Kernel, typename Approximate_kernel::Construct_bbox_2, typename Exact_kernel::Construct_bbox_2>             Construct_bbox_2;
  typedef Lazy_construction_bbox<Kernel, typename Approximate_kernel::Construct_bbox_3, typename Exact_kernel::Construct_bbox_3>             Construct_bbox_3;
  typedef Lazy_cartesian_const_iterator_2<Kernel, typename Approximate_kernel::Construct_cartesian_const_iterator_2, typename Exact_kernel::Construct_cartesian_const_iterator_2>   Construct_cartesian_const_iterator_2;
  typedef Lazy_cartesian_const_iterator_3<Kernel, typename Approximate_kernel::Construct_cartesian_const_iterator_3, typename Exact_kernel::Construct_cartesian_const_iterator_3>   Construct_cartesian_const_iterator_3;

  typedef CGAL::CartesianKernelFunctors::Compute_approximate_squared_length_3<Kernel>  Compute_approximate_squared_length_3;
  typedef CGAL::CartesianKernelFunctors::Compute_approximate_area_3<Kernel>  Compute_approximate_area_3;

  // typedef void Compute_z_3; // to detect where .z() is called
  // typedef void Construct_point_3; // to detect where the ctor is called

  Assign_2
  assign_2_object() const
  { return Assign_2(); }

  Assign_3
  assign_3_object() const
  { return Assign_3(); }

  Construct_bbox_2
  construct_bbox_2_object() const
  { return Construct_bbox_2(); }

  Construct_bbox_3
  construct_bbox_3_object() const
  { return Construct_bbox_3(); }

  Construct_cartesian_const_iterator_2
  construct_cartesian_const_iterator_2_object() const
  { return Construct_cartesian_const_iterator_2(); }

  Construct_cartesian_const_iterator_3
  construct_cartesian_const_iterator_3_object() const
  { return Construct_cartesian_const_iterator_3(); }

  Compute_approximate_squared_length_3
  compute_approximate_squared_length_3_object() const
  { return Compute_approximate_squared_length_3(); }

  Compute_approximate_area_3
  compute_approximate_area_3_object() const
  { return Compute_approximate_area_3(); }
}; // end class Lazy_kernel_base<EK_, AK_, E2A_, Kernel_2>

#ifndef CGAL_LAZY_KERNEL_USE_STATIC_FILTERS_BY_DEFAULT
#  ifdef CGAL_NO_STATIC_FILTERS
#    define CGAL_LAZY_KERNEL_USE_STATIC_FILTERS_BY_DEFAULT false
#  else 
#    define CGAL_LAZY_KERNEL_USE_STATIC_FILTERS_BY_DEFAULT true
#  endif
#endif

template <class Exact_kernel, class Approximate_kernel, class E2A>
struct Lazy_kernel_without_type_equality
  : public Lazy_kernel_base< Exact_kernel, Approximate_kernel, E2A, Lazy_kernel_without_type_equality<Exact_kernel,Approximate_kernel, E2A> >
{};

template <class Exact_kernel,
	  class Approximate_kernel = Simple_cartesian<Interval_nt_advanced>,
          class E2A = Cartesian_converter<Exact_kernel, Approximate_kernel>,
          bool UseStaticFilters = CGAL_LAZY_KERNEL_USE_STATIC_FILTERS_BY_DEFAULT >
struct Lazy_kernel
  : public Type_equality_wrapper<
             Lazy_kernel_base< Exact_kernel, Approximate_kernel, E2A,
                               Lazy_kernel<Exact_kernel, Approximate_kernel, E2A, UseStaticFilters> >,
             Lazy_kernel<Exact_kernel, Approximate_kernel, E2A, UseStaticFilters> >
{};

template <class Exact_kernel, class Approximate_kernel, class E2A>
struct Lazy_kernel<Exact_kernel, Approximate_kernel, E2A, true>
  : public internal::Static_filters<
      Type_equality_wrapper<
        Lazy_kernel_base< Exact_kernel, Approximate_kernel, E2A, Lazy_kernel<Exact_kernel, Approximate_kernel, E2A, true> > ,
        Lazy_kernel<Exact_kernel, Approximate_kernel, E2A, true> >, false >
{
// WARNING: If you change the definition of Lazy_kernel, then you need to
// change also the definition of Epeck in
// <CGAL/Exact_predicate_exact_constructions_kernel.h>.
};

} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_LAZY_KERNEL_H

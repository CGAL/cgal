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

#include <boost/mpl/if.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>


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

template<bool has_result_type, typename F>
struct Maybe_result_type {
  // This is incredibly evil. It relies on the fact that we always
  // have a type in the primary result template but it looks like the
  // only way out given the current design.
  typedef typename F::template result<void>::type type;
};

template<typename F>
struct Maybe_result_type<true, F> {
  typedef typename F::result_type type;
};

// goes through the standard process of selecting the right
// Lazy_something after the kind of the return type has been
// determined
template<typename T, typename AK, typename EK, typename Kernel, typename AKC, typename EKC>
struct Standard_pick {
  typedef typename boost::remove_cv< typename boost::remove_reference< typename T::type >::type >::type T_;
  typedef typename boost::mpl::if_< boost::is_same< T_, typename AK::FT  >,
                                    Lazy_construction_nt<Kernel, AKC, EKC>,
                                    typename boost::mpl::if_< boost::is_same< T_, Object >,
                                                              Lazy_construction_object<Kernel, AKC, EKC>,
                                                              Lazy_construction<Kernel, AKC, EKC> >::type
  >::type type;
};
} // internal

// Exact_kernel = exact kernel that will be made lazy
// Kernel = lazy kernel

// the Generic base simplies applies the generic magic functor stupidly.
// then the real base fixes up a few special cases.
template < typename EK_, typename AK_, typename E2A_, typename Kernel >
class Lazy_kernel_generic_base
  // : public Filtered_kernel_base<EK_>
    // TODO : Static_filters_base too ?  Check performance
{
public:

  typedef AK_   Approximate_kernel;
  typedef EK_   Exact_kernel;
  typedef E2A_  E2A;

  // 3 synonyms identical to Filtered_kernel (TODO : cleanup !)
  typedef AK_   FK;
  //typedef E2A_  C2F;

  typedef Approx_converter<Kernel, Approximate_kernel>   C2F;
  typedef Exact_converter<Kernel, Exact_kernel>    C2E;
  // Note: Approx_converter and Exact_converter are defined in <CGAL/Lazy.h>

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

  // What to do with the tag ?
  // Probably this should not exist, should it ?
  // struct filter_tag{};
  // typedef filter_tag                                     Kernel_tag;
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


  // We don't touch the predicates.
  // FIXME TODO : better use a layer of Filtered_kernel on top of everything,
  //              so that semi-static filters are used as well (?).
#define CGAL_Kernel_pred(P, Pf)  \
    typedef Filtered_predicate<typename Exact_kernel::P, typename Approximate_kernel::P, C2E, C2F> P; \
    P Pf() const { return P(); }

    // We change the constructions.
#ifdef CGAL_INTERSECT_WITH_ITERATORS_2
#define CGAL_Kernel_cons(C, Cf) \
    typedef typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C, typename Approximate_kernel::Intersect_with_iterators_2>, \
                                     Lazy_intersect_with_iterators<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                     typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Bbox_2>, \
                                                              Lazy_construction_bbox<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                                              typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, typename Approximate_kernel::FT>,\
                                                                                       Lazy_construction_nt<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>,\
                                                                                       typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Object >,\
                                                                                                                Lazy_construction_object<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>,\
                                                                                                                Lazy_construction<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C> >::type >::type > ::type > ::type C; \
    C Cf() const { return C(); }

  CGAL_Kernel_cons(Intersect_with_iterators_2,
		   intersect_with_iterators_2_object)
#else
#define CGAL_Kernel_cons(C, Cf) \
  typedef typename internal::Standard_pick< internal::Maybe_result_type< internal::Has_result_type< typename Approximate_kernel::C >::value, \
                                                                         typename Approximate_kernel::C >, \
                                            Approximate_kernel, Exact_kernel, Kernel, typename Approximate_kernel::C, typename Exact_kernel::C \
                                            >::type C;                  \
  C Cf() const { return C(); }

#endif //CGAL_INTERSECT_WITH_ITERATORS_2

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

#if 0
    // We change the constructions.
#ifdef CGAL_INTERSECT_WITH_ITERATORS_2
#define CGAL_Kernel_cons(C, Cf) \
    typedef typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C, typename Approximate_kernel::Intersect_with_iterators_2>, \
                                     Lazy_intersect_with_iterators<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                     typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Bbox_2>, \
                                                              Lazy_construction_bbox<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                                              typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, typename Approximate_kernel::FT>,\
                                                                                       Lazy_construction_nt<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>,\
                                                                                       typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Object >,\
                                                                                                                Lazy_construction_object<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>,\
                                                                                                                Lazy_construction<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C> >::type >::type > ::type > ::type C; \
    C Cf() const { return C(); }

  CGAL_Kernel_cons(Intersect_with_iterators_2,
		   intersect_with_iterators_2_object)
#else
#define CGAL_Kernel_cons(C, Cf) \
    typedef typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C, typename Approximate_kernel::Construct_cartesian_const_iterator_2>, \
                                     Lazy_cartesian_const_iterator_2<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                     typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C, typename Approximate_kernel::Construct_cartesian_const_iterator_3>, \
                                                              Lazy_cartesian_const_iterator_3<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                                              typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Bbox_2>, \
                                                                                       Lazy_construction_bbox<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                                                                       typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Bbox_3>, \
                                                                                                                Lazy_construction_bbox<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>, \
                                                                                                                typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, typename Approximate_kernel::FT>,\
                                                                                                                                         Lazy_construction_nt<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>,\
                                                                                                                                         typename boost::mpl::if_<boost::is_same<typename Approximate_kernel::C::result_type, Object >,\
                                                                                                                                                                  Lazy_construction_object<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C>,\
                                                                                                                                                                  Lazy_construction<Kernel,typename Approximate_kernel::C, typename Exact_kernel::C> >::type >::type >::type > ::type > ::type > ::type C; \
    C Cf() const { return C(); }

#endif //CGAL_INTERSECT_WITH_ITERATORS_2

#endif // 0

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

#endif // CGAL_LAZY_KERNEL_H

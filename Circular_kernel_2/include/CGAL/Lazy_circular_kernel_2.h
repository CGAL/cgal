// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Andreas Fabri, Constantinos Tsirogiannis

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)


#ifndef CGAL_LAZY_CIRCULAR_KERNEL_H
#define CGAL_LAZY_CIRCULAR_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Circular_kernel_type_equality_wrapper.h>
#include <CGAL/Algebraic_kernel_converter.h>
#include <CGAL/Circular_kernel_converter.h>
#include <CGAL/Lazy.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Lazy_circular_kernel_constructions.h>
#include <CGAL/Root_of_traits.h>
#include <boost/mpl/if.hpp>
#include <CGAL/Kernel/Type_mapper.h>


namespace CGAL {

// Some specializations (move to another file ?)
template < typename K1, typename K2 >
struct Type_mapper < typename K1::Circular_arc_2, K1, K2 >
{
  typedef typename K2::Circular_arc_2 type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Circular_arc_point_2, K1, K2 >
{
  typedef typename K2::Circular_arc_point_2 type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Line_arc_2, K1, K2 >
{
  typedef typename K2::Line_arc_2 type;
};


// EK = exact kernel that will be made lazy
// Kernel = lazy kernel
template < typename EK_, typename AK_, typename E2A_, typename Kernel >
class Lazy_circular_kernel_base
//  : public EK::template Base<Kernel>::Type
#if 1
  : public Lazy_kernel_base< EK_, AK_, E2A_,
                             Lazy_kernel<EK_, AK_, E2A_> >
#endif
{
  //    typedef typename EK::template Base<Kernel>::Type   Kernel_base;
    // Hardcoded for now.
  //    typedef Simple_cartesian<Interval_nt_advanced>   AK; // An optimiser
public:
  typedef AK_ AK;
  typedef typename EK_::Algebraic_kernel  Algebraic_kernel;
  //  typedef Simple_cartesian<Interval_nt<> >   AK;
  typedef EK_   EK;
  typedef E2A_  E2A;

  typedef Approx_converter<Kernel, AK>   C2F;
  typedef Exact_converter<Kernel, EK>    C2E;

#if 0
  typedef typename EK::Algebraic_kernel Al_EK;
  typedef typename AK::Algebraic_kernel Al_AK;
  typedef Algebraic_kernel_converter<Al_EK,Al_AK,
                                     To_interval<typename EK::RT>, 
				     To_interval<typename EK::Root_of_2> > Al_K_converter;
	
  typedef typename EK::Linear_kernel Ln_EK;
  typedef typename AK::Linear_kernel Ln_AK;
  typedef Cartesian_converter<Ln_EK, Ln_AK,
                                To_interval<typename Ln_EK::RT> > Ln_converter;
				
  // typedef Circular_kernel_converter<EK,AK,Ln_converter,Al_K_converter> E2A;
				
  //typedef Lazy_kernel<Ln_EK,Ln_AK>   Linear_kernel;
   typedef Lazy_kernel<Ln_EK,Ln_AK,E2A>   Linear_kernel0;
   typedef typename Linear_kernel0::template Base<Kernel>::Type Linear_kernel; // should not be needed.

#endif

   template < typename Kernel2 >
   struct Base { typedef Lazy_circular_kernel_base<EK, AK, E2A, Kernel2>  Type; };

    // What to do with the tag ?
    // Probably this should not exist, should it ?
    // struct filter_tag{};
    // typedef filter_tag                                     Kernel_tag;
    // typedef typename CK::Kernel_tag                       Kernel_tag;
    // typedef typename CK::Rep_tag                          Rep_tag;


    //Please remove this if you consider it to be sloppy
      struct  Lazy_tag{};
      typedef Lazy_tag                     Definition_tag;
      typedef typename EK::Definition_tag  Circular_tag;
    ////////////////////////////////////////////////////

    // Types
  typedef CGAL::Lazy_exact_nt<typename EK::FT>  FT;
  typedef FT RT;
  typedef typename Root_of_traits<RT>::RootOf_2  Root_of_2;

//  typedef Lazy<typename AK::Algebraic_kernel::Root_of_2, typename EK::Algebraic_kernel::Root_of_2, 
//            typename EK::FT, Al_K_converter> Root_of_2;

  typedef FT Cartesian_coordinate_type;
  typedef RT Homogeneous_coordinate_type; //Nothing to do with the circular_k
  					    //just for the time being

  typedef CGAL::Object Object_2;
  typedef CGAL::Object Object_3;

  typedef Lazy<typename AK::Circular_arc_2, typename EK::Circular_arc_2, 
  	       typename EK::FT, E2A>                                        Circular_arc_2;
	       
  typedef Lazy<typename AK::Circular_arc_point_2, 
               typename EK::Circular_arc_point_2, 
	       typename EK::FT, E2A>                                        Circular_arc_point_2;
	       
  typedef Lazy<typename AK::Line_arc_2, typename EK::Line_arc_2, 
  	       typename EK::FT, E2A>                                        Line_arc_2;

  //  typedef Root_for_circles_2_2<RT> Root_for_circles_2_2;
  typedef int Root_for_circles_2_2; // should AK be filtered ?

    // We don't touch the predicates.
#define CGAL_Circular_Kernel_pred(P, Pf)  \
    typedef Filtered_predicate< typename EK::P, typename AK::P, C2E, C2F> P; \
    P Pf() const { return P(); }


    // We change the constructions.

#define CGAL_Circular_Kernel_cons(C, Cf) \
    typedef typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Intersect_2>, \
                                     Lazy_construct_intersections_2 <Kernel,typename AK::C, typename EK::C>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Make_x_monotone_2>, \
                                     Lazy_make_x_monotone_2 <Kernel,typename AK::C, typename    EK::C>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Make_xy_monotone_2>, \
                                     Lazy_make_x_monotone_2 <Kernel,typename AK::C, typename    EK::C>, \
	    typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Split_2>, \
                                     Lazy_functor_2_2 <Kernel, typename AK::C, typename EK::C>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C::result_type, Bbox_2>, \
                                     Lazy_construction_bbox<Kernel,typename AK::C, typename EK::C>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C::result_type, typename AK::RT>,\
                                     Lazy_construction_nt<Kernel,typename AK::C, typename EK::C>,\
            typename boost::mpl::if_<boost::is_same<typename AK::C::result_type, Object >,\
                                     Lazy_construction_object<Kernel,typename AK::C, typename EK::C>,\
         Lazy_construction<Kernel,typename AK::C, typename EK::C> >::type >::type >::type >::type >::type >::type >::type C; \
    C Cf() const { return C(); }



#include <CGAL/Circular_kernel_2/interface_macros.h>

};

template <class EK, class AK, class E2A>
struct Lazy_circular_kernel_adaptor
  : public Lazy_circular_kernel_base< EK, AK, E2A, Lazy_circular_kernel_adaptor<EK, AK, E2A> >
{};

template <class EK, class AK, class E2A>
struct Lazy_circular_kernel_without_type_equality
  : public Lazy_circular_kernel_base< EK, AK, E2A, Lazy_circular_kernel_without_type_equality<EK, AK, E2A> >
{};

template <class EK, 
  class AK = void ,//Circular_kernel<Cartesian<Interval_nt_advanced>,
                   //        Algebraic_kernel_2_2<Interval_nt_advanced> >,
  class E2A = Circular_kernel_converter<EK,AK,
                 Cartesian_converter< EK, AK, To_interval<typename EK::RT> >,
                 Algebraic_kernel_converter< typename EK::Algebraic_kernel, typename AK::Algebraic_kernel,
                 To_interval<typename EK::RT>, To_interval<typename EK::Root_of_2> > > >
struct Lazy_circular_kernel_2
  : public Circular_kernel_type_equality_wrapper< 
             Lazy_circular_kernel_base< EK, AK, E2A, Lazy_circular_kernel_2<EK, AK, E2A> >,
             Lazy_circular_kernel_2<EK, AK, E2A> >
{};

} //namespace CGAL

#endif // CGAL_LAZY_CIRCULAR_KERNEL_H

// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Sylvain Pion, Constantinos Tsirogiannis

#ifndef CGAL_LAZY_CURVED_KERNEL_H
#define CGAL_LAZY_CURVED_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Curved_kernel_type_equality_wrapper.h>
#include <CGAL/Algebraic_kernel_converter.h>
#include <CGAL/Curved_kernel_converter.h>
#include <CGAL/Lazy.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Lazy_curved_kernel_constructions.h>
#include <CGAL/Root_of_2.h>
#include <boost/mpl/if.hpp>


CGAL_BEGIN_NAMESPACE


// EK = exact kernel that will be made lazy
// Kernel = lazy kernel
template < typename EK_, typename AK_, typename Kernel >
class Lazy_curved_kernel_base
//  : public EK::template Base<Kernel>::Type
{
  //    typedef typename EK::template Base<Kernel>::Type   Kernel_base;
    // Hardcoded for now.
  //    typedef Simple_cartesian<Interval_nt_advanced>   AK; // An optimiser
public:
  typedef AK_ AK;
  //  typedef Simple_cartesian<Interval_nt<> >   AK;
  typedef EK_   EK;
  typedef typename EK::Algebraic_kernel Al_EK;
  typedef typename AK::Algebraic_kernel Al_AK;
  typedef Algebraic_kernel_converter<Al_EK,Al_AK,
                                     To_interval<typename EK::RT>, 
				     To_interval<typename EK::Root_of_2> > Al_K_converter;
	
  typedef typename EK::Linear_kernel Ln_EK;
  typedef typename AK::Linear_kernel Ln_AK;
  typedef Cartesian_converter<Ln_EK, Ln_AK,
                                To_interval<typename Ln_EK::RT> > Ln_converter;
				
  typedef Curved_kernel_converter<EK,AK,Ln_converter,Al_K_converter> E2A;
				
   typedef Lazy_kernel<Ln_EK,Ln_AK>   Linear_kernel;


    template < typename Kernel2 >
    struct Base { typedef Lazy_curved_kernel_base<EK, AK, Kernel2>  Type; };

    // What to do with the tag ?
    // Probably this should not exist, should it ?
    // struct filter_tag{};
    // typedef filter_tag                                     Kernel_tag;
    // typedef typename CK::Kernel_tag                       Kernel_tag;
    // typedef typename CK::Rep_tag                          Rep_tag;

    // Types
  typedef CGAL::Lazy_exact_nt<typename EK::FT>  FT;
  typedef FT RT;
  typedef FT Cartesian_coordinate_type;
  typedef RT Homogeneous_coordinate_type; //Nothing to do with the curved_k
  					    //just for the time being

  typedef CGAL::Object Object_2;
  typedef CGAL::Object Object_3;

  typedef Lazy<typename Al_AK::Root_of_2, typename Al_EK::Root_of_2, 
               typename EK::FT, Al_K_converter> Root_of_2;
 
  typedef typename Lazy_kernel<Ln_EK,Ln_AK>::Point_2  Point_2;
  typedef typename Lazy_kernel<Ln_EK,Ln_AK>::Segment_2  Segment_2;
  typedef typename Lazy_kernel<Ln_EK,Ln_AK>::Line_2   Line_2;
  typedef typename Lazy_kernel<Ln_EK,Ln_AK>::Circle_2 Circle_2;
  typedef typename Lazy_kernel<Ln_EK,Ln_AK>::Conic_2  Conic_2;
  typedef typename Lazy_kernel<Ln_EK,Ln_AK>::Cartesian_const_iterator_2  
					     Cartesian_const_iterator_2;

  
  typedef Lazy<typename AK::Circular_arc_2, typename EK::Circular_arc_2, 
  	       typename EK::FT, E2A>                                        Circular_arc_2;
	       
  typedef Lazy<typename AK::Circular_arc_point_2, 
               typename EK::Circular_arc_point_2, 
	       typename EK::FT, E2A>                                        Circular_arc_point_2;
	       
  typedef Lazy<typename AK::Line_arc_2, typename EK::Line_arc_2, 
  	       typename EK::FT, E2A>                                        Line_arc_2;

	       
  

    // We don't touch the predicates.
#define CGAL_Curved_Kernel_pred(P, Pf)  \
    typedef Filtered_predicate< typename EK::P, typename AK::P,Exact_converter, Approx_converter> P; \
    P Pf() const { return P(); }


    // We change the constructions.

#define CGAL_Curved_Kernel_cons(C, Cf) \
    typedef typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Intersect_2>, \
                                     Lazy_construct_intersections_2 <Kernel,AK,EK,typename AK::C, typename EK::C, \
				     typename EK::FT, E2A>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Make_x_monotone_2>, \
                                     Lazy_make_x_monotone_2 <Kernel,AK,EK,typename AK::C, typename    EK::C, \
				     typename EK::FT, E2A>, \
	    typename boost::mpl::if_<boost::is_same<typename AK::C, typename AK::Split_2>, \
                                     Lazy_functor_2_2 <Kernel,AK,EK,typename AK::C, typename EK::C, \
				     typename EK::FT, E2A>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C::result_type, Bbox_2>, \
                                     Lazy_construction_bbox<AK,EK,typename AK::C, typename EK::C, typename EK::FT, E2A>, \
            typename boost::mpl::if_<boost::is_same<typename AK::C::result_type, Interval_nt<true> >,\
                                     Lazy_construction_nt<AK,EK,typename AK::C, typename EK::C, typename EK::FT, E2A>,\
            typename boost::mpl::if_<boost::is_same<typename AK::C::result_type, Object >,\
                                     Lazy_construction_object<Kernel,AK,EK,typename AK::C, typename EK::C, typename EK::FT, E2A>,\
                                     Lazy_construction<AK,EK,typename AK::C, typename EK::C, typename EK::FT, E2A> >::type >::type > ::type > ::type > ::type > ::type C; \
    C Cf() const { return C(); }



#include <CGAL/Curved_kernel/interface_macros.h>

};

template <class EK, class AK>
struct Lazy_curved_kernel_adaptor
  : public Lazy_curved_kernel_base< EK, AK, Lazy_curved_kernel_adaptor<EK,AK> >
{};

template <class EK, class AK>
struct Lazy_curved_kernel_without_type_equality
  : public Lazy_curved_kernel_base< EK, AK, Lazy_curved_kernel_without_type_equality<EK,AK> >
{};

template <class EK, class AK = void>
struct Lazy_curved_kernel
  : public Curved_kernel_type_equality_wrapper< 
             Lazy_curved_kernel_base< EK, AK, Lazy_curved_kernel<EK, AK> >,
             Lazy_curved_kernel<EK, AK> >
{};

CGAL_END_NAMESPACE

#endif // CGAL_LAZY_CURVED_KERNEL_H






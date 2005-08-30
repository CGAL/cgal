// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Circular_kernel.h

#ifndef CGAL_CIRCULAR_KERNEL_H
#define CGAL_CIRCULAR_KERNEL_H

#include <CGAL/Curved_kernel/Circular_arc_endpoint_2.h>
#include <CGAL/Curved_kernel/Circular_arc_2.h>
#include <CGAL/Curved_kernel/Line_arc_2.h>
#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_endpoint_2.h>
#include <CGAL/Line_arc_2.h>

#include <CGAL/Curved_kernel/function_objects_on_circle_2.h>
#include <CGAL/global_functions_on_circle_2.h>

#include <CGAL/Curved_kernel/function_objects_polynomial_circular.h>
#include <CGAL/global_functions_on_circular_arcs_2.h>

namespace CGAL {
namespace CGALi {

template < class CurvedKernel, class LinearKernel >
struct Curved_kernel_base: public LinearKernel
// takes classes in internal sub-namespace
{
  typedef CGALi::Circular_arc_2<CurvedKernel>               Circular_arc_2;
  typedef CGALi::Circular_arc_point_2<CurvedKernel>      Circular_arc_point_2;
  typedef CGALi::Line_arc_2<CurvedKernel>                   Line_arc_2;
  
  
   
  
  #define CGAL_Curved_Kernel_pred(Y,Z) typedef CircularFunctors::Y<CurvedKernel> Y; \
                              Y Z() const { return Y(); }
  #define CGAL_Curved_Kernel_cons(Y,Z) CGAL_Curved_Kernel_pred(Y,Z)

  #include <CGAL/Curved_kernel/interface_macros.h>
};

} // namespace CGALi

template < class LinearKernel, class AlgebraicKernel >
struct Curved_kernel
  :  // there should be a derivation from
    // LinearKernel::Kernel_base<Self> to have types equalities for
    // the Linearkernel types
    public CGALi::Curved_kernel_base<Curved_kernel<LinearKernel,AlgebraicKernel>,LinearKernel >
{
  typedef LinearKernel                                    Linear_kernel;
  typedef AlgebraicKernel                                 Algebraic_kernel;

  typedef Curved_kernel<LinearKernel,AlgebraicKernel>      Self;
  typedef CGALi::Curved_kernel_base<Self,LinearKernel>     Kernel_base;

  typedef typename LinearKernel::RT                       RT;
  typedef typename LinearKernel::FT                       FT;

  typedef typename Algebraic_kernel::Root_of_2            Root_of_2;
  typedef typename Algebraic_kernel::Root_for_circles_2_2  Root_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_for_circles_2_2
                                                             Polynomial_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_1_2
                                                  Polynomial_1_2;

  // public classes
  typedef CGAL::Object Object_2;
  typedef CGAL::Object Object_3;
  
  typedef typename Linear_kernel::Line_2                  Line_2;
  typedef typename Linear_kernel::Circle_2                Circle_2;
  typedef typename Linear_kernel::Conic_2                 Conic_2;
  typedef typename Linear_kernel::Point_2                 Point_2;


  typedef CGAL::Circular_arc_2<Self>                      Circular_arc_2;
  typedef CGAL::Circular_arc_point_2<Self>                Circular_arc_point_2;
  typedef CGAL::Line_arc_2<Self>                          Line_arc_2;




};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_H

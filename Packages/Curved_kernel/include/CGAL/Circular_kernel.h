// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

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

#include <CGAL/Curved_kernel_type_equality_wrapper.h>

namespace CGAL {
namespace CGALi {

template < class CurvedKernel, class LinearKernelBase >
struct Curved_kernel_base: public LinearKernelBase
{
  typedef CGALi::Circular_arc_2<CurvedKernel>         Circular_arc_2;
  typedef CGALi::Circular_arc_point_2<CurvedKernel>   Circular_arc_point_2;
  typedef CGALi::Line_arc_2<CurvedKernel>             Line_arc_2;
  
  #define CGAL_Curved_Kernel_pred(Y,Z) \
    typedef CircularFunctors::Y<CurvedKernel> Y; \
    Y Z() const { return Y(); }
  #define CGAL_Curved_Kernel_cons(Y,Z) CGAL_Curved_Kernel_pred(Y,Z)

  #include <CGAL/Curved_kernel/interface_macros.h>
};

} // namespace CGALi

template < class LinearKernel, class AlgebraicKernel >
struct Curved_kernel
  : public Curved_kernel_type_equality_wrapper
  <
  CGALi::Curved_kernel_base
  < Curved_kernel<LinearKernel,AlgebraicKernel>,
    typename LinearKernel::template 
    Base<Curved_kernel<LinearKernel,AlgebraicKernel> >::Type 
  >,
  Curved_kernel<LinearKernel,AlgebraicKernel> 
  >
{
  typedef Curved_kernel<LinearKernel,AlgebraicKernel>      Self;

  typedef typename LinearKernel::template 
  Base<Curved_kernel<LinearKernel,AlgebraicKernel> >::Type Linear_kernel;
  typedef AlgebraicKernel                                  Algebraic_kernel;

  // for Lazy hexagons/bbox kernels
  // Please remove this if you consider it to be sloppy
  struct Curved_tag{};
  typedef Curved_tag Definition_tag;
  //

  typedef typename LinearKernel::RT                       RT;
  typedef typename LinearKernel::FT                       FT;

  typedef typename Algebraic_kernel::Root_of_2            Root_of_2;
  typedef typename Algebraic_kernel::Root_for_circles_2_2 Root_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_for_circles_2_2
                                                    Polynomial_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_1_2 Polynomial_1_2;

//   typedef CGAL::Object Object_2;
//   typedef CGAL::Object Object_3;
  
};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_H

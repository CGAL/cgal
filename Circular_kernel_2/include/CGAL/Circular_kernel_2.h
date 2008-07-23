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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_2_H
#define CGAL_CIRCULAR_KERNEL_2_H

#include <CGAL/Circular_kernel_2/Circular_arc_point_2.h>
#include <CGAL/Circular_kernel_2/Circular_arc_2.h>
#include <CGAL/Circular_kernel_2/Line_arc_2.h>
#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_point_2.h>
#include <CGAL/Line_arc_2.h>

#include <CGAL/Circular_kernel_2/function_objects_on_circle_2.h>
#include <CGAL/global_functions_on_circle_2.h>

#include <CGAL/Circular_kernel_2/function_objects_polynomial_circular.h>
#include <CGAL/global_functions_on_circular_arcs_2.h>

#include <CGAL/Circular_kernel_type_equality_wrapper.h>

namespace CGAL {
namespace CGALi {

template < class CircularKernel, class LinearKernelBase >
struct Circular_kernel_base_ref_count: public LinearKernelBase
{
  typedef CGALi::Circular_arc_2<CircularKernel>         Circular_arc_2;
  typedef CGALi::Circular_arc_point_2<CircularKernel>   Circular_arc_point_2;
  typedef CGALi::Line_arc_2<CircularKernel>             Line_arc_2;

  // The mecanism that allows to specify reference-counting or not.
  template < typename T >
  struct Handle { typedef Handle_for<T>    type; };
  
  #define CGAL_Circular_Kernel_pred(Y,Z) \
    typedef CircularFunctors::Y<CircularKernel> Y; \
    Y Z() const { return Y(); }
  #define CGAL_Circular_Kernel_cons(Y,Z) CGAL_Circular_Kernel_pred(Y,Z)

  #include <CGAL/Circular_kernel_2/interface_macros.h>
};

} // namespace CGALi

template < class LinearKernel, class AlgebraicKernel >
struct Circular_kernel_2
  : public Circular_kernel_type_equality_wrapper
  <
  CGALi::Circular_kernel_base_ref_count
  < Circular_kernel_2<LinearKernel,AlgebraicKernel>,
    typename LinearKernel:: template 
    Base<Circular_kernel_2<LinearKernel,AlgebraicKernel> >::Type 
  >,
  Circular_kernel_2<LinearKernel,AlgebraicKernel> 
  >
{
  typedef Circular_kernel_2<LinearKernel,AlgebraicKernel>      Self;

  typedef typename LinearKernel::template 
  Base<Circular_kernel_2<LinearKernel,AlgebraicKernel> >::Type Linear_kernel;
  typedef AlgebraicKernel                                  Algebraic_kernel;

  // for Lazy hexagons/bbox kernels
  // Please remove this if you consider it to be sloppy
  struct Circular_tag{};
  typedef Circular_tag Definition_tag;
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

#endif // CGAL_CIRCULAR_KERNEL_2_H

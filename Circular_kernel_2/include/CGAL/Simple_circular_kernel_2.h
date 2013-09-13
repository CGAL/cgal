// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SIMPLE_CIRCULAR_KERNEL_2_H
#define CGAL_SIMPLE_CIRCULAR_KERNEL_2_H

#if 0

#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_point_2.h>
#include <CGAL/Line_arc_2.h>
#include <CGAL/Circular_kernel_2/Circular_arc_point_2.h>
#include <CGAL/Circular_kernel_2/Circular_arc_2.h>
#include <CGAL/Circular_kernel_2/Line_arc_2.h>

#include <CGAL/Circular_kernel_type_equality_wrapper.h>

#include <CGAL/Circular_kernel_2/function_objects_polynomial_circular.h>
#include <CGAL/global_functions_circular_kernel_2.h>

#include <CGAL/Circular_kernel_2/function_objects_on_line_2.h>

#include <CGAL/Circular_kernel_2/function_objects_on_circle_2.h>

namespace CGAL {

namespace internal {

template < class CircularKernel, class LinearKernelBase, class AlgebraicKernel >
struct Circular_kernel_base_ref_count: public LinearKernelBase
{
  typedef internal::Circular_arc_2<CircularKernel>         Circular_arc_2;
  typedef internal::Circular_arc_point_2<CircularKernel>   Circular_arc_point_2;
  typedef internal::Line_arc_2<CircularKernel>             Line_arc_2;
  typedef LinearKernelBase                              Linear_kernel;
  typedef AlgebraicKernel                               Algebraic_kernel;
  typedef typename Algebraic_kernel::Root_of_2            Root_of_2;
  typedef typename Algebraic_kernel::Root_for_circles_2_2 Root_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_for_circles_2_2
                                                    Polynomial_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_1_2 Polynomial_1_2;
  typedef typename Linear_kernel::RT                       RT;
  typedef typename Linear_kernel::FT                       FT;

  // The mechanism that allows to specify reference-counting or not.
  template < typename T >
  struct Handle { typedef T    type; };

  template < typename Kernel2 >
  struct Base { typedef Circular_kernel_base_ref_count<Kernel2, LinearKernelBase, AlgebraicKernel>  Type; };  

  #define CGAL_Circular_Kernel_pred(Y,Z) \
    typedef CircularFunctors::Y<CircularKernel> Y; \
    Y Z() const { return Y(); }
  #define CGAL_Circular_Kernel_cons(Y,Z) CGAL_Circular_Kernel_pred(Y,Z)

  #include <CGAL/Circular_kernel_2/interface_macros.h>
};

} // namespace internal

template < class LinearKernel, class AlgebraicKernel >
struct Circular_kernel_2
  : public Circular_kernel_type_equality_wrapper
  <
  internal::Circular_kernel_base_ref_count
  < Circular_kernel_2<LinearKernel,AlgebraicKernel>,
    typename LinearKernel:: template 
    Base<Circular_kernel_2<LinearKernel,AlgebraicKernel> >::Type,
    AlgebraicKernel 
  >,
  Circular_kernel_2<LinearKernel,AlgebraicKernel>
  >
{
  // for Lazy hexagons/bbox kernels
  // Please remove this if you consider it to be sloppy
  struct Circular_tag{};
  typedef Circular_tag Definition_tag;
  //  
};

} //namespace CGAL

#endif

#endif // CGAL_SIMPLE_CIRCULAR_KERNEL_2_H

// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_2_H
#define CGAL_CIRCULAR_KERNEL_2_H

#include <CGAL/license/Circular_kernel_2.h>


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
  typedef internal::Circular_arc_2_base<CircularKernel>         Circular_arc_2;
  typedef internal::Circular_arc_point_2_base<CircularKernel>   Circular_arc_point_2;
  typedef internal::Line_arc_2_base<CircularKernel>             Line_arc_2;
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
  struct Handle { typedef Handle_for<T>    type; };

  template < typename Kernel2 >
  struct Base {
    typedef typename LinearKernelBase::template Base<Kernel2>::Type ReboundLK;
    typedef Circular_kernel_base_ref_count<Kernel2,
                                           ReboundLK,
                                           AlgebraicKernel>  Type;
  };

  #define CGAL_Circular_Kernel_pred(Y,Z) \
    typedef CircularFunctors::Y<CircularKernel> Y; \
    Y Z() const { return Y(); }
  #define CGAL_Circular_Kernel_cons(Y,Z) CGAL_Circular_Kernel_pred(Y,Z)

  #include <CGAL/Circular_kernel_2/interface_macros.h>

  typedef LinearFunctors::Construct_line_2<CircularKernel> Construct_line_2;
  Construct_line_2 construct_line_2_object() const { return Construct_line_2(); }
};

} // namespace internal

template < class LinearKernel, class AlgebraicKernel >
struct Circular_kernel_2
  : public Circular_kernel_type_equality_wrapper
     < internal::Circular_kernel_base_ref_count
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

#endif // CGAL_CIRCULAR_KERNEL_2_H

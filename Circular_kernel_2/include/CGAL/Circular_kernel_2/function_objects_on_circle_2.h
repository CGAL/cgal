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


#ifndef CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H
#define CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/Circular_kernel_2/internal_functions_on_circle_2.h>

#include <CGAL/Circular_kernel_2/function_objects_on_line_2.h>
// to be removed when CGAL::Kernel has a Get_equation

namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  class Construct_circle_2
    // : public CK::Linear_kernel::Construct_circle_2
  {
    typedef typename CK::FT                         FT;
    typedef typename CK::Circle_2                   Circle_2;
    typedef typename CK::Circular_arc_2             Circular_arc_2;

    typedef typename CK::Linear_kernel::Construct_circle_2 Linear_Construct_circle_2;

  public:
    // using Linear_Construct_circle_2::operator();

    template <class... Args>
    decltype(auto)
    operator()(const Args&... args) const
    { return Linear_Construct_circle_2()(args...); }

    Circle_2
    operator() ( const typename CK::Polynomial_for_circles_2_2 &eq ) {
      return construct_circle_2<CK>(eq);
    }

    decltype(auto)
    operator() (const Circular_arc_2 & a) const {
      return (a.rep().supporting_circle());
    }

  };

} // namespace CircularFunctors

} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H

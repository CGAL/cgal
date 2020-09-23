// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado,
//             Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_3_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {
  namespace SphericalFunctors {

  // we can optimize those functions by comparing
  // the references before doing the comparison
  // as in CK
  template < class SK >
  inline
  Comparison_result
  compare_x(const typename SK::Circular_arc_point_3 &p0,
            const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   Algebraic_kernel;
    return Algebraic_kernel().compare_x_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result
  compare_y(const typename SK::Circular_arc_point_3 &p0,
            const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   Algebraic_kernel;
    return Algebraic_kernel().compare_y_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result
  compare_z(const typename SK::Circular_arc_point_3 &p0,
            const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   Algebraic_kernel;
    return Algebraic_kernel().compare_z_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result
  compare_xy(const typename SK::Circular_arc_point_3 &p0,
             const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   Algebraic_kernel;
    return Algebraic_kernel().compare_xy_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result
  compare_xyz(const typename SK::Circular_arc_point_3 &p0,
              const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   Algebraic_kernel;
    return Algebraic_kernel().compare_xyz_object()(p0.coordinates(), p1.coordinates());
  }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_3_H

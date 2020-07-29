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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_PLANE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_PLANE_3_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {
  namespace SphericalFunctors {

    template < class SK >
    typename SK::Plane_3
    construct_plane_3(const typename SK::Polynomial_1_3 &eq)
    {
      typedef typename SK::Plane_3 Plane_3;
      return Plane_3(eq.a(),eq.b(),eq.c(),eq.d());
    }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_PLANE_3_H

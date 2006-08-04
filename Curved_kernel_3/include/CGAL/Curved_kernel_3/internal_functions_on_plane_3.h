// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_PLANE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_PLANE_3_H

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

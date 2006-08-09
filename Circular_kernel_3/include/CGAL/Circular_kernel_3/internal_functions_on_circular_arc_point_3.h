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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_POINT_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_POINT_3_H

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
      equal( const typename SK::Circular_arc_point_3 &p1,
             const typename SK::Circular_arc_point_3 &p2)
    {
      return p1.rep() == p2.rep();
    }
        
    
  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_POINT_3_H

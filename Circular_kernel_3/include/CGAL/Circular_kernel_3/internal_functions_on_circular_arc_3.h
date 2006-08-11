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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
    equal( const typename SK::Circular_arc_3 &c1,
           const typename SK::Circular_arc_3 &c2)
    {
      return c1.rep() == c2.rep();
    }

    template< class SK>
    inline
    Sign
    compute_sign_of_cross_product(const typename SK::Root_of_2 &x1, 
                                  const typename SK::Root_of_2 &y1,
                                  const typename SK::Root_of_2 &z1,
                                  const typename SK::Root_of_2 &x2, 
                                  const typename SK::Root_of_2 &y2,
                                  const typename SK::Root_of_2 &z2) {
      typedef typename SK::Root_of_2 Root_of_2;
      const Root_of_2 cx = y1 * z2 - z1 * y2;
      const Root_of_2 cy = z1 * x2 - x1 * z2;
      const Root_of_2 cz = x1 * y2 - y1 * x2;
      if(!is_zero(cx)) return sign(cx);
      if(!is_zero(cy)) return sign(cy);
      return sign(cz);
    }

    template< class SK>
    inline
    Sign
    compute_sign_of_cross_product(const typename SK::FT &x1, 
                                  const typename SK::FT &y1,
                                  const typename SK::FT &z1,
                                  const typename SK::FT &x2, 
                                  const typename SK::FT &y2,
                                  const typename SK::FT &z2) {
      typedef typename SK::FT FT;
      const FT cx = y1 * z2 - z1 * y2;
      const FT cy = z1 * x2 - x1 * z2;
      const FT cz = x1 * y2 - y1 * x2;
      if(!is_zero(cx)) return sign(cx);
      if(!is_zero(cy)) return sign(cy);
      return sign(cz);
    }

    template< class SK>
    inline
    Sign
    compute_sign_of_cross_product(const typename SK::Circular_arc_point_3 &p1, 
                                  const typename SK::Circular_arc_point_3 &p2,
                                  const typename SK::Point_3 &c) {
      return compute_sign_of_cross_product<SK>(p1.x()-c.x(),
                                               p1.y()-c.y(),
                                               p1.z()-c.z(),
                                               p2.x()-c.x(),
                                               p2.y()-c.y(),
                                               p2.z()-c.z());
    }

    template< class SK>
    inline
    Sign
    compute_sign_of_cross_product(const typename SK::Point_3 &p1, 
                                  const typename SK::Point_3 &p2,
                                  const typename SK::Point_3 &c) {
      return compute_sign_of_cross_product<SK>(p1.x()-c.x(),
                                               p1.y()-c.y(),
                                               p1.z()-c.z(),
                                               p2.x()-c.x(),
                                               p2.y()-c.y(),
                                               p2.z()-c.z());
    }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H


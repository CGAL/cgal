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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_3_H

namespace CGAL {
  namespace SphericalFunctors {

    template < class SK >
    typename SK::Line_3
    construct_line_3(const typename SK::Polynomials_for_line_3 &eq)
    {
      typedef typename SK::Line_3 Line_3;
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Vector_3 Vector_3;
      return Line_3(Point_3(eq.b1(),eq.b2(),eq.b3()),
                    Vector_3(eq.a1(),eq.a2(),eq.a3()));
    }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_3_H

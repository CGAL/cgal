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

#include <CGAL/Circular_kernel_3/internal_function_has_on_spherical_kernel.h>

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
    equal( const typename SK::Circular_arc_3 &c1,
           const typename SK::Circular_arc_3 &c2)
    {
      return c1.rep() == c2.rep();
    }

    template <class SK>
    inline
    bool
    do_overlap(const typename SK::Circular_arc_3 &c1,
               const typename SK::Circular_arc_3 &c2,
               const bool known_equal_supporting_circle = false)
    { 
      if(!known_equal_supporting_circle) {
        if(!non_oriented_equal<SK>(c1.supporting_circle(), 
                                   c2.supporting_circle()))
          return false;
      }
      if(c1.rep().is_full()) return true;
      if(c2.rep().is_full()) return true;
      if((SK().has_on_3_object()(c1,c2.target(),true)) || 
         (SK().has_on_3_object()(c1,c2.source(),true))) return true;
      return SK().has_on_3_object()(c2,c1.source(),true);
    }

    template < class SK >
    void
    split(const typename SK::Circular_arc_3 &c,
	  const typename SK::Circular_arc_point_3 &p,
	  typename SK::Circular_arc_3 &c1,
	  typename SK::Circular_arc_3 &c2)
    {
      // The point must be on the circular arc 
      CGAL_kernel_precondition(SK().has_on_3_object()(c, p));
      typedef typename SK::Circular_arc_3  Circular_arc_3;
      // It doesn't make sense to split an arc on an extremity
      CGAL_kernel_precondition(c.source() != p);
      CGAL_kernel_precondition(c.target() != p);
      const Circular_arc_3 &rc1 = 
        Circular_arc_3(c.supporting_circle(), c.source(), p);
      const Circular_arc_3 &rc2 = 
        Circular_arc_3(c.supporting_circle(), p, c.target());
      if ( SK().compare_xyz_3_object()(rc1.source(), rc2.source()) != 
           SMALLER) {
        c1 = rc2; c2 = rc1;
      } else { c1 = rc1; c2 = rc2; }
    }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H


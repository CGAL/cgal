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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_3_H

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
    typedef typename SK::Algebraic_kernel   AK;
    return AK().compare_x_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result 
  compare_y(const typename SK::Circular_arc_point_3 &p0,
            const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   AK;
    return AK().compare_y_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result 
  compare_z(const typename SK::Circular_arc_point_3 &p0,
            const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   AK;
    return AK().compare_z_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result 
  compare_xy(const typename SK::Circular_arc_point_3 &p0,
             const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   AK;
    return AK().compare_xy_object()(p0.coordinates(), p1.coordinates());
  }

  template < class SK >
  inline
  Comparison_result 
  compare_xyz(const typename SK::Circular_arc_point_3 &p0,
              const typename SK::Circular_arc_point_3 &p1)
  {
    typedef typename SK::Algebraic_kernel   AK;
    return AK().compare_xyz_object()(p0.coordinates(), p1.coordinates());
  }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_3_H

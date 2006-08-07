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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
    equal( const typename SK::Line_arc_3 &l1,
           const typename SK::Line_arc_3 &l2)
    {
      return l1.rep() == l2.rep();
    }

    template <class SK>
    inline
    bool
    do_overlap(const typename SK::Line_arc_3 &l1,
               const typename SK::Line_arc_3 &l2)
    { 
      if (!non_oriented_equal<SK>(l1.supporting_line(),
                                  l2.supporting_line())) 
        return false;

      return SK().compare_xyz_3_object()(l1.higher_xyz_extremity(), 
                             l2.lower_xyz_extremity()) >= 0
          && SK().compare_xyz_3_object()(l1.lower_xyz_extremity(), 
                             l2.higher_xyz_extremity()) <= 0;
    }

    template < class SK >
    void
    split(const typename SK::Line_arc_3 &l,
	  const typename SK::Circular_arc_point_3 &p,
	  typename SK::Line_arc_3 &l1,
	  typename SK::Line_arc_3 &l2)
    {
      typedef typename SK::Line_arc_3  Line_arc_3;
      CGAL_kernel_precondition(SK().has_on_3_object()(l, p));
      // It doesn't make sense to split an arc on an extremity
      CGAL_kernel_precondition(l.source() != p);
      CGAL_kernel_precondition(l.target() != p);
      if(SK().compare_xyz_3_object()(l.source(),p) == SMALLER) {
        l1 = Line_arc_3(l.supporting_line(),l.source(),p);
        l2 = Line_arc_3(l.supporting_line(),p,l.target());
      } else {
        l1 = Line_arc_3(l.supporting_line(),p,l.target());
        l2 = Line_arc_3(l.supporting_line(),l.source(),p);
      } 
    }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_ARC_3_H

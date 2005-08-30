// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Curved_kernel/Circular_arc_endpoint_2.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_POINT_2_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_POINT_2_H

#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Cartesian.h>
#include <iostream>
#include <cassert>
#include <CGAL/Curved_kernel/Debug_id.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Interval_arithmetic.h>

#include <CGAL/global_functions_on_circle_2.h>
//#include <CGAL/global_functions_on_roots_and_polynomials_2_2.h> 
// fixme, devrait
// appeler fonction de global_functions_on_circular_arcs

namespace CGAL {
namespace CGALi {

  template <class CK >
  class Circular_arc_point_2
    : public Debug_id<>
  {
    typedef typename CK::FT                      FT;
    typedef typename CK::Root_of_2               Root_of_2;

  public: // fixme ?

    typedef typename CK::Root_for_circles_2_2 Root_for_circles_2_2;
    
    Circular_arc_point_2() {}

    Circular_arc_point_2(const Root_for_circles_2_2 & np)
      :  _p(np){}

    const Root_of_2 & x() const { return _p.x(); }
    const Root_of_2 & y() const { return _p.y(); }

    CGAL::Bbox_2 bbox() const
    {
       std::pair<double,double> ix=to_interval(x()),
                                iy=to_interval(y());

        return CGAL::Bbox_2(ix.first,iy.first,
                      ix.second,iy.second);
    }

    const Root_for_circles_2_2 & coordinates() const { return _p; }

  private:
    Root_for_circles_2_2 _p;
  };

  template < class CK >
  std::ostream&
  operator<<(std::ostream &os, const Circular_arc_point_2<CK> &p)
  {
    return os << p.x()[0] << " " 
	      << p.x()[1] << " " 
	      << p.x()[2] << " " 
	      << p.x().is_smaller() << " "
	      << p.y()[0] << " " 
	      << p.y()[1] << " " 
	      << p.y()[2] << " "
	      << p.y().is_smaller() << " "; 
  }
  
  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_point_2<CK> &p)
  {
    typedef typename CK::Root_of_2               Root_of_2;
    typedef typename CK::Root_for_circles_2_2 Root_for_circles_2_2;

    typename Root_of_2::RT x1, x2, x3;
    typename Root_of_2::RT y1, y2, y3;
    bool b1, b2;
    is >> x1 >> x2 >> x3 >> b1 >> y1 >> y2 >> y3 >> b2 ;
    if (is)
      p = Circular_arc_point_2<CK>(Root_for_circles_2_2(Root_of_2(x3, x2, x1, b1),
						      Root_of_2(y3, y2, y1, b2)));
    return is;
  }

 template < typename CK >
  std::ostream &
  print(std::ostream & os, const Circular_arc_point_2<CK> &p)
  {
   return os << "CirclArcEndPoint_2(" << p.id() << std::endl
	      << p.x() << ", " << p.y() << ')';
  }

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_ENDPOINT_2_H

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

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_ENDPOINT_2_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_ENDPOINT_2_H

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
  class Circular_arc_endpoint_2
    : public Debug_id<>
  {
    typedef typename CK::Circle_2                Circle_2;
    typedef typename CK::Line_2                  Line_2;
    typedef typename CK::FT                      FT;
    typedef typename CK::Root_of_2               Root_of_2;

  public: // fixme ?
    typedef typename CGAL::Simple_cartesian<Root_of_2>::Point_2
                                                 Numeric_point_2;

#if 0
    static const Circular_arc_endpoint_2 def; // default constructed instance.

    struct pipo {}; // Used to select the private constructor.

    Circular_arc_endpoint_2(pipo)
    {
      typedef typename CK::Linear_kernel::Point_2       Pt_2;
      *this = Circular_arc_endpoint_2(Circle_2(Pt_2(0, 1), 1),
				      Circle_2(Pt_2(1, 0), 1), true);
    }
#endif

    // An auxiliary function
    // (we could abstract it away to a multivariate solver in a further step)
    static
    Numeric_point_2
    intersect(const Circle_2 & c1, const Circle_2 & c2, const bool b)
    {

      typedef std::vector<CGAL::Object >
	solutions_container;

      solutions_container solutions;
      CGAL::construct_intersections_2<CK>
	( c1, c2, std::back_inserter(solutions) );
      typename solutions_container::iterator it = solutions.begin();

      assert( it != solutions.end() ); // the circles intersect
      std::pair<typename CK::Circular_arc_endpoint_2, uint> result;
      std::cout << assign(result, *it) << std::endl;
      if ( result.second == 2 ) // double solution
	return result.first._p;
      if (b) {
	return result.first._p;
      }
      ++it;
      bool tmp = assign(result, *it);
      assert(tmp);
      return result.first._p;

 //       FT dx = c2.center().x() - c1.center().x();
//       FT dy = c2.center().y() - c1.center().y();
//       FT px = c2.center().x() + c1.center().x();
//       FT py = c2.center().y() + c1.center().y();
//       FT dx2 = CGAL::square(dx);
//       FT dy2 = CGAL::square(dy);
//       FT dist2 = dx2 + dy2;
//       FT rx1 = c1.squared_radius() - CGAL::square(c1.center().x());
//       FT ry1 = c1.squared_radius() - CGAL::square(c1.center().y());
//       FT rx2 = c2.squared_radius() - CGAL::square(c2.center().x());
//       FT ry2 = c2.squared_radius() - CGAL::square(c2.center().y());
//       
//       FT drx = rx2 - rx1;
//       FT dry = ry2 - ry1;

//       bool low_y = ( CGAL::sign(dx) *
//                      CGAL::sign(dy) <= 0 ) // !!! case = 0 ???
// 	           ? b : (! b );

//       return Numeric_point_2(
//              make_root_of_2(4*dist2,
//                             4*(dx*drx-px*dy2),
//                             CGAL::square(drx) - dy2*(2*(rx1+rx2)-dy2),
//                             b),
//              make_root_of_2(4*dist2,
//                             4*(dy*dry-py*dx2),
//                             CGAL::square(dry) - dx2*(2*(ry1+ry2)-dx2),
//                             low_y));

    }

  public:

    // FIXME : this is one more source of inefficiency...
    // It seems we need a valid default constructed endpoint,
    // Otherwise point_reflect_x_and_y is not happy (because the 2 circles
    // don't intersect necessarily like they should)...
    // So let's construct the origin this way.
    // We copy the static instance to avoid too many useless computations.
    Circular_arc_endpoint_2()
      // : _c0(def._c0), _c1(def._c1), _b(def._b), _p(def._p) {}
    {
      // let's try with a static instance here instead.
      typedef typename CK::Linear_kernel::Point_2       Pt_2;
      static const Circular_arc_endpoint_2 def = 
	Circular_arc_endpoint_2(Circle_2(Pt_2(0, 1), 1),
				Circle_2(Pt_2(1, 0), 1), true, Numeric_point_2());
      *this = def;
    }

    Circular_arc_endpoint_2(const Circle_2 & c1, const Circle_2 & c2, 
			    const bool b)
      : _c0(c1), _c1(c2), _b(b), _p(intersect(c1, c2, b)) {}

    Circular_arc_endpoint_2(const Circle_2 & c1, const Circle_2 & c2,
			    const bool b, const Numeric_point_2 & np)
      : _c0(c1), _c1(c2), _b(b), _p(np){}

    const Root_of_2 & x() const { return _p.x(); }
    const Root_of_2 & y() const { return _p.y(); }

    const Circle_2 & circle(int i) const { return i==0 ? _c0 : _c1; }

    bool is_left() const { return _b; }

    Bbox_2 bbox() const
    {
       std::pair<double,double> ix=to_interval(x()),
                                iy=to_interval(y());

       return Bbox_2 (ix.first,iy.first,
                      ix.second,iy.second);
    }

  private:
    Circle_2  _c0, _c1;
    bool      _b;
    Numeric_point_2 _p;
  };

  // template < class CK >
  // const Circular_arc_endpoint_2<CK> Circular_arc_endpoint_2<CK>::def = Circular_arc_endpoint_2<CK>(pipo());

  template < class CK >
  std::ostream&
  operator<<(std::ostream &os, const Circular_arc_endpoint_2<CK> &p)
  {
    return os << "CirclArcEndPoint_2(" << p.id() << std::endl
              << p.circle(0) << " ," << std::endl
              << p.circle(1) << " ," << std::endl
              << p.is_left() << " , " << p.x() << ", " << p.y() << ')';
  }

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_ENDPOINT_2_H

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

// file : include/CGAL/Curved_kernel/Circular_arc_2.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_2_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_2_H

#include <CGAL/Curved_kernel/Debug_id.h>
#include <CGAL/global_functions_on_circular_arcs_2.h>
#include <CGAL/intersections.h>

namespace CGAL {
namespace CGALi {

  template <class CK >
  class Circular_arc_2
    : public Debug_id<>
  {
    typedef typename CK::FT                        FT;
    typedef typename CK::RT                        RT;
    typedef typename CK::Linear_kernel::Point_2    Point_2;
    typedef typename CK::Line_2                    Line_2;
    typedef typename CK::Circle_2                  Circle_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;
    typedef typename CK::Root_of_2                 Root_of_2;

  public:

    Circular_arc_2() {}

    // Full circle
    Circular_arc_2(const Circle_2 &c)
      : _swap_source(false)
    {
       // Define a circle intersecting c in the 2 vertical tangent
       // points.
       Circle_2 c1 (Point_2(c.center().x(), c.center().y()-1),
                    c.squared_radius()+1);

       _begin = Circular_arc_endpoint_2(c, c1, true);
       _end   = Circular_arc_endpoint_2(c, c1, true);
    }

    Circular_arc_2(const Circle_2 &support,
                   const Line_2 &l1, bool b1,
                   const Line_2 &l2, bool b2)
      : _swap_source(false)
    {
      Point_2 center1 (support.center().x() + l1.a()/2,
                       support.center().y() + l1.b()/2);

      FT sqr1 = support.squared_radius() + l1.c()
                - CGAL::square(support.center().x())
                - CGAL::square(support.center().y())
                + CGAL::square(center1.x())
                + CGAL::square(center1.y());

      Circle_2 c1 (center1, sqr1);

      Point_2 center2 (support.center().x() + l2.a()/2,
                       support.center().y() + l2.b()/2);

      FT sqr2 = support.squared_radius() + l2.c()
                - CGAL::square(support.center().x())
                - CGAL::square(support.center().y())
                + CGAL::square(center2.x())
                + CGAL::square(center2.y());

      Circle_2 c2 (center2, sqr2);

      *this = Circular_arc_2(support, c1, b1, c2, b2);

      assert(do_intersect(support, c1));
      assert(do_intersect(support, c2));
    }

    Circular_arc_2(const Circle_2 &c, 
		   const Circle_2 &c1, const bool b_1,
		   const Circle_2 &c2, const bool b_2)
      : _begin(c, c1, b_1), _end(c, c2, b_2), _swap_source(false) {}


    // constructs a circular arc that is the arc included in A
    // having same (b) endpoint as A (true == _begin, false == _end)
    // but whose (!b) endpoint is the intersection of A with ccut given 
    // by b_cut
    Circular_arc_2(const Circular_arc_2 &A, const bool b,
		   const Circle_2 &ccut, const bool b_cut)
      : _swap_source(A._swap_source)
    {
        assert(A.is_x_monotone());
        assert(do_intersect(A.supporting_circle(), ccut));

        Circular_arc_endpoint_2 new_p (A.supporting_circle(), ccut, b_cut);

        // assert(point_in_range(A, new_p));
        assert(A.on_upper_part() ==
              (CGAL::compare(new_p.y(), A.center().y()) >= 0));

        if (b) {
          _begin  = A._begin;
          _end    = new_p;
        }
        else {
          _begin  = new_p;
          _end    = A._end;
        }
      }

    // Constructs an arc supported by Circle_2(begin, middle, end),
    // with _begin == begin, _end == end.
    // (middle is not necessarily on the arc)
    Circular_arc_2(const Point_2 &begin,
                   const Point_2 &middle,
                   const Point_2 &end)
      : _swap_source(false)
    {
      assert(!collinear(begin, middle, end));

      Circle_2 c  (begin, middle, end);
      Line_2   l1 (begin, middle);
      Line_2   l2 (middle, end);
      *this = Circular_arc_2(c, l1, compare_xy(begin, middle) < 0,
                                l2, compare_xy(end,   middle) < 0);
    }
    
    // Constructs an arc supported by Circle_2,
    // with _begin == begin, _end == end.
    Circular_arc_2(const Circle_2 &support,
                   const Point_2 &source,
                   const Point_2 &target)
      : _swap_source(false)
    {
      
      assert(support.has_on_boundary(source));
      assert(support.has_on_boundary(target));
        
      Line_2   l (source, target);
      bool source_is_left =  (compare_xy(source, target) < 0) ;
      *this = Circular_arc_2(support, l, source_is_left ,l, !source_is_left);
    }

  private:

    // The arc goes from _begin to _end in the positive order
    // If _begin == _end, then it's the full circle
    Circular_arc_endpoint_2  _begin, _end;
    // If swap_source is false, then source == _begin and target == _end.
    // Otherwise it's the reverse.
    bool        _swap_source;

  public:

    bool swap_source() const
    {
      return _swap_source;
    }

    void swap_source_and_target()
    {
      _swap_source = !_swap_source;
    }

    const Circular_arc_endpoint_2 & source() const
    {
      return _swap_source ? _end : _begin;
    }

    const Circular_arc_endpoint_2 & target() const
    {
      return _swap_source ? _begin : _end;
    }

    const Circular_arc_endpoint_2 & left() const
    {
      assert(is_x_monotone());
      assert(on_upper_part() ? compare_xy(_end,_begin)<0
	     : compare_xy(_begin,_end)<0);
      return on_upper_part() ? _end : _begin;
    }

    const Circular_arc_endpoint_2 & right() const
    {
      assert(is_x_monotone());
      assert(on_upper_part() ? compare_xy(_end,_begin)<0
	     : compare_xy(_begin,_end)<0);
      return on_upper_part() ? _begin : _end;
    }

    const Circular_arc_endpoint_2 & begin() const
    {
      return _begin;
    }

    const Circular_arc_endpoint_2 & end() const
    {
      return _end;
    }

    bool is_x_monotone() const
    {
      int cmp_begin = CGAL::compare(_begin.y(), center().y());
      int cmp_end   = CGAL::compare(_end.y(),   center().y());

      // XXX : be careful, this may be surprising if the return value
      // is not -1/1 but some random int...
      if (cmp_begin == opposite(cmp_end) && cmp_begin != 0)
        return false;

      int cmp_x = compare_x(_begin, _end);

      // Is the arc on the upper part ?
      if (cmp_begin > 0 || cmp_end > 0)
        return cmp_x > 0;

      // Is the arc on the lower part ?
      if (cmp_begin < 0 || cmp_end < 0)
        return cmp_x < 0;

      // There remains the case :
      assert(cmp_begin == 0 && cmp_end == 0);

      return cmp_x != 0; // full circle or half circle.
    }

    bool on_upper_part() const
      // check whether the endpoints are above or below the center
    { 
      assert(is_x_monotone());

      int begin_y = CGAL::compare(_begin.y(), supporting_circle().center().y());
      int end_y   = CGAL::compare(_end.y(),   supporting_circle().center().y());

      if (begin_y == 0 && end_y == 0)
        return compare_x(_begin, _end) > 0;

      return begin_y > 0 || end_y > 0;
    }

    const Circle_2 & supporting_circle() const           
    {
       assert(_begin.circle(0) == _end.circle(0));
       return _begin.circle(0);
    }

    const Point_2 & center() const           
    {
       return supporting_circle().center();
    }

    const FT & squared_radius() const           
    {
       return supporting_circle().squared_radius();
    }

    // Until curves_compare_y_at_x() is fixed in PM.
    double approximate_y_at(const Circular_arc_endpoint_2 &p) const
    {
      assert(is_x_monotone());
      double x = CGAL::to_double(p.x());
      double tmp = std::sqrt(CGAL::to_double(supporting_circle()
                                                .squared_radius())
                   - CGAL::square(x -
                         CGAL::to_double(supporting_circle().center().x())));
      if (on_upper_part())
        return CGAL::to_double(supporting_circle().center().y()) + tmp;
      else
        return CGAL::to_double(supporting_circle().center().y()) - tmp;
    }
  };

  template < typename CK >
  std::ostream &
  operator<<(std::ostream & os, const Circular_arc_2<CK> &a)
  {
    // The output format is :
    // - supporting circle
    // - circle c1
    // - bool b1
    // - circle c2
    // - bool b2
    return os << a.supporting_circle() << " "
	      << a.begin().circle(1) << " " << a.begin().is_left() << " "
	      << a.end().circle(1)   << " " << a.end().is_left() << " ";
  }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_2<CK> &a)
  {
    typename CK::Circle_2 s, c1, c2;
    bool b1, b2;
    is >> s >> c1 >> b1 >> c2 >> b2;
    if (is)
      a = Circular_arc_2<CK>(s, c1, b1, c2, b2);
    return is;
  }

  template < typename CK >
  std::ostream &
  print(std::ostream & os, const Circular_arc_2<CK> &a)
  {
    return os << "Circular_arc_2( " << a.id() << std::endl
              << "left : " << a.left() << " , " << std::endl
              << "right : " << a.right() << " , " << std::endl
	      << "upper part : " << a.on_upper_part() << std::endl
              << "  [[ approximate circle is (x,y,r) : "
              << CGAL::to_double(a.supporting_circle().center().x()) << "  "
              << CGAL::to_double(a.supporting_circle().center().y()) << "  "
              << std::sqrt(CGAL::to_double(a.supporting_circle()
                                            .squared_radius()))
              << " ]]" << std::endl;
  }

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_2_H

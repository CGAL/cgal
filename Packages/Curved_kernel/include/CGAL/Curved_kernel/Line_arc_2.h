#ifndef CGAL_CURVED_KERNEL_LINE_ARC_2_H
#define CGAL_CURVED_KERNEL_LINE_ARC_2_H

#include <CGAL/Curved_kernel/Debug_id.h>
#include <CGAL/global_functions_on_line_2.h>
#include <CGAL/global_functions_on_circle_2.h>
#include <CGAL/global_functions_on_line_arcs_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Algebraic_kernel/internal_functions_on_roots_and_polynomial_1_2_and_2_2.h>
#include <CGAL/Curved_kernel/internal_functions_on_line_2.h>
#include <CGAL/Curved_kernel/internal_functions_on_line_arc_2.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {
namespace CGALi {

  template <class CK >
  class Line_arc_2
    : public Debug_id<>
  {
    typedef typename CK::FT                        FT;
    typedef typename CK::RT                        RT;
    typedef typename CK::Linear_kernel::Point_2    Point_2;
    typedef typename CK::Line_2                    Line_2;
    typedef typename CK::Circle_2                  Circle_2;
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Circular_arc_point_2   Circular_arc_point_2;
    typedef typename CK::Root_of_2                 Root_of_2;
    typedef typename CK::Segment_2                 Segment_2;

  public:
    //typedef typename CGAL::Simple_cartesian<Root_of_2>::Point_2
    //                                             Numeric_point_2;
    typedef typename CK::Algebraic_kernel::Root_for_circles_2_2 
      Root_for_circles_2_2;
    
    static
      Circular_arc_point_2
      intersect(const Line_2 & l, const Circle_2 & c, const bool b)
    {
      
      typedef std::vector<CGAL::Object >
	solutions_container;
      
      solutions_container solutions;
      CGAL::LinearFunctors::intersect_2<CK>
	( l, c, std::back_inserter(solutions) );
      typename solutions_container::iterator it = solutions.begin();
      
      CGAL_kernel_precondition( it != solutions.end() ); 
      // the circles intersect
      
      const std::pair<typename CK::Circular_arc_point_2, uint> *result;
      result = CGAL::object_cast< 
        std::pair<typename CK::Circular_arc_point_2, uint> >(&(*it));
      if ( result->second == 2 ) // double solution
	return result->first;
      if (b) return result->first;
      ++it;
      result = CGAL::object_cast< 
        std::pair<typename CK::Circular_arc_point_2, uint> >(&(*it));
      return result->first;
    }


  public:
    Line_arc_2() {}
     
    Line_arc_2(const Line_2 &support,
	       const Circle_2 &c1,const bool b1,
	       const Circle_2 &c2,const bool b2)
      :_support(support)
    {
      _begin = intersect(support, c1, b1);
      _end = intersect(support, c2, b2);


    }


    Line_arc_2(const Line_2 &support,
	       const Line_2 &l1,
	       const Line_2 &l2)
      :_support(support)
    {
      CGAL_kernel_precondition(do_intersect(support, l1));
      CGAL_kernel_precondition(do_intersect(support, l2));
      typedef typename Root_of_2::RT RT_2;
      //Voir pour mettre une assertion au assign
      Object obj = intersection(support, l1);
      const Point_2 *pt = CGAL::object_cast<Point_2>(&obj);
      _begin = Circular_arc_point_2(Root_for_circles_2_2(Root_of_2(pt->x()),
						       Root_of_2(pt->y())));
      obj = intersection(support, l2);
      const Point_2 *pt2 = CGAL::object_cast<Point_2>(&obj);
      _end = Circular_arc_point_2(Root_for_circles_2_2(Root_of_2(pt2->x()),
						     Root_of_2(pt2->y())));
    }
    
    Line_arc_2(const Line_2 &support,
	       const Circular_arc_point_2 &p1,
	       const Circular_arc_point_2 &p2)
      :_support(support)
    {
      //Verifier si p1 et p2 sont sur la line
      _begin = p1;
      _end = p2;
    }

    Line_arc_2(const Segment_2 &s)
      :_support(s.supporting_line())
      ,_begin(s.source())
      ,_end(s.target())
    {}
    

    Line_arc_2(const Point_2 &p1,
	       const Point_2 &p2)
    {
      _support = Line_2(p1, p2);
      _begin = Circular_arc_point_2(Root_for_circles_2_2(Root_of_2(p1.x()),
						       Root_of_2(p1.y())));
      _end = Circular_arc_point_2(Root_for_circles_2_2(Root_of_2(p2.x()),
						     Root_of_2(p2.y())));
    }

  private:
    
    Line_2 _support;
    Circular_arc_point_2 _begin, _end;


  public :
    
    const Line_2 & supporting_line() const
    {
      return _support;
    }
    
    const Circular_arc_point_2 & left() const
    {
      return compare_xy(_begin, _end) < 0 ? _begin : _end;
    }

    const Circular_arc_point_2 & right() const
    {
      return compare_xy(_begin, _end) < 0 ? _end : _begin;
    }

    const Circular_arc_point_2 & source() const
    {
      return _begin;
    }

    const Circular_arc_point_2 & target() const
    {
      return _end;
    } 
    
    CGAL::Bbox_2 bbox() const
    {
      return _begin.bbox() + _end.bbox();
    }

  };
  
  template < typename CK >
    std::ostream &
    operator<<(std::ostream & os, const Line_arc_2<CK> &a)
    {
      
      return os << a.supporting_line() << " "
		<< a.source() << " "
		<< a.target() << " ";
    }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Line_arc_2<CK> &a)
  {
    typename CK::Line_2 l;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> l >> p1 >> p2 ;
    if (is)
      a = Line_arc_2<CK>(l, p1, p2);
    return is;
  }


 } // namespace CGALi
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_LINE_ARC_2_H

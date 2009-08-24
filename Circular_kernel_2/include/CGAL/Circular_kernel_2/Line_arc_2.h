// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_LINE_ARC_2_H
#define CGAL_CIRCULAR_KERNEL_LINE_ARC_2_H

#include <CGAL/global_functions_circular_kernel_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Algebraic_kernel_for_circles/internal_functions_on_roots_and_polynomial_1_2_and_2_2.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_line_2.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_line_arc_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Circular_kernel_2/Circular_arc_2.h>

namespace CGAL {
namespace internal {

  template <class CK >
  class Line_arc_2
  {
    typedef typename CK::FT                        FT;
    typedef typename CK::RT                        RT;
    typedef typename CK::Point_2                   Point_2;
    typedef typename CK::Line_2                    Line_2;
    typedef typename CK::Circle_2                  Circle_2;
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Circular_arc_point_2      Circular_arc_point_2;
    typedef typename CK::Root_of_2                 Root_of_2;
    typedef typename CK::Segment_2                 Segment_2;
    typedef struct bit_field {
      unsigned char begin_less_xy_than_end:2;
    } bit_field;

    private:
    // set flags to 0
    // when 1 bit -> 0 = false, 1 = true
    // when 2 bits -> 0 = don_know, 1 = false
    //                              2 = true
    void reset_flags() const {
      flags.begin_less_xy_than_end = 0;
    }

  public:
    //typedef typename CGAL::Simple_cartesian<Root_of_2>::Point_2
    //                                             Numeric_point_2;
    typedef typename CK::Root_for_circles_2_2 
      Root_for_circles_2_2;
    
    static
      Circular_arc_point_2
      intersect(const Line_2 & l, const Circle_2 & c, const bool b)
    {
      
      typedef std::vector<CGAL::Object >
	solutions_container;
      
      solutions_container solutions;
      CGAL::CircularFunctors::intersect_2<CK>
	( l, c, std::back_inserter(solutions) );
      typename solutions_container::iterator it = solutions.begin();
      
      CGAL_kernel_precondition( it != solutions.end() ); 
      // the circles intersect
      
      const std::pair<typename CK::Circular_arc_point_2, unsigned> *result;
      result = CGAL::object_cast< 
        std::pair<typename CK::Circular_arc_point_2, unsigned> >(&(*it));
      if ( result->second == 2 ) // double solution
	return result->first;
      if (b) return result->first;
      ++it;
      result = CGAL::object_cast< 
        std::pair<typename CK::Circular_arc_point_2, unsigned> >(&(*it));
      return result->first;
    }


  public:
    Line_arc_2() 
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    : id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id()) 
#endif 
    {}
     
    Line_arc_2(const Line_2 &support,
	       const Circle_2 &c1,const bool b1,
	       const Circle_2 &c2,const bool b2)
      :_support(support)
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    ,id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id()) 
#endif 
    {
      _begin = intersect(support, c1, b1);
      _end = intersect(support, c2, b2);
      reset_flags();
    }


    Line_arc_2(const Line_2 &support,
	       const Line_2 &l1,
	       const Line_2 &l2)
      :_support(support)
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    ,id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id()) 
#endif 
    {
      CGAL_kernel_precondition(do_intersect(support, l1));
      CGAL_kernel_precondition(do_intersect(support, l2));
      //typedef typename Root_of_2::RT RT_2;
      //Voir pour mettre une assertion au assign
      Object obj = intersection(support, l1);
      const Point_2 *pt = CGAL::object_cast<Point_2>(&obj);
      _begin = Circular_arc_point_2(*pt);
      obj = intersection(support, l2);
      const Point_2 *pt2 = CGAL::object_cast<Point_2>(&obj);
      _end = Circular_arc_point_2(*pt2);
      reset_flags();
    }
    
    Line_arc_2(const Line_2 &support,
	       const Circular_arc_point_2 &p1,
	       const Circular_arc_point_2 &p2)
      :_support(support)
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    ,id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id()) 
#endif 
    {
      //Verifier si p1 et p2 sont sur la line
      _begin = p1;
      _end = p2;
      reset_flags();
    }

    Line_arc_2(const Segment_2 &s)
      :_support(s.supporting_line())
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    ,id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id()) 
#endif 
    {
      _begin = Circular_arc_point_2(s.source());
      _end = Circular_arc_point_2(s.target());
      reset_flags();
    }
    

    Line_arc_2(const Point_2 &p1,
	       const Point_2 &p2)
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    : id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id()) 
#endif 
    {
      _support = Line_2(p1, p2);
      _begin = Circular_arc_point_2(p1);
      _end = Circular_arc_point_2(p2);
      reset_flags();
    }

  private:
    
    Line_2 _support;
    Circular_arc_point_2 _begin, _end;
    mutable bit_field flags;

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES 
    mutable unsigned int id_of_my_supporting_line;

    unsigned int line_number() const {
      return id_of_my_supporting_line;
    }

    void set_line_number(unsigned int i) const {
      id_of_my_supporting_line = i;
    }
#endif

  private: //(some useful functions)

    bool begin_less_xy_than_end() const {
      if(flags.begin_less_xy_than_end == 0) {
        if(compare_xy(_begin, _end) < 0)
          flags.begin_less_xy_than_end = 2;
        else flags.begin_less_xy_than_end = 1;
      } return flags.begin_less_xy_than_end == 2;
    }

  public :

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES 
    template < class T >
    static bool find_intersection_circle_line(
      const Circular_arc_2& c, 
      const Line_arc_2& l, 
      T& res) {
      if(c.id_of_my_supporting_circle == 0) return false;
      return Circular_arc_2::circle_table.template find<T>(c.id_of_my_supporting_circle, 
                                  l.id_of_my_supporting_line, 
                                  res);
    }

    template < class T >
    static void put_intersection_circle_line(const Circular_arc_2& c, 
      const Line_arc_2& l,
      const T& res) {
      Circular_arc_2::circle_table.template put<T>(c.circle_number(), 
                          l.line_number(), 
                          res);
    }
#endif
    
    const Line_2 & supporting_line() const
    {
      return _support;
    }
    
    const Circular_arc_point_2 & left() const
    {
      return begin_less_xy_than_end() ? _begin : _end;
    }

    const Circular_arc_point_2 & right() const
    {
      return begin_less_xy_than_end() ? _end : _begin;
    }

    const Circular_arc_point_2 & source() const
    {
      return _begin;
    }

    const Circular_arc_point_2 & target() const
    {
      return _end;
    } 
    
    bool is_vertical() const
    {
      return supporting_line().is_vertical();
    }
    
    CGAL::Bbox_2 bbox() const
    {
      return _begin.bbox() + _end.bbox();
    }

  };
  
  /* template < typename CK > */
/*     std::ostream & */
/*     operator<<(std::ostream & os, const Line_arc_2<CK> &a) */
/*     { */
      
/*       return os << a.supporting_line() << " " */
/* 		<< a.source() << " " */
/* 		<< a.target() << " "; */
/*     } */

/*   template < typename CK > */
/*   std::istream & */
/*   operator>>(std::istream & is, Line_arc_2<CK> &a) */
/*   { */
/*     typename CK::Line_2 l; */
/*     typename CK::Circular_arc_point_2 p1; */
/*     typename CK::Circular_arc_point_2 p2; */
/*     is >> l >> p1 >> p2 ; */
/*     if (is) */
/*       a = Line_arc_2<CK>(l, p1, p2); */
/*     return is; */
/*   } */


 } // namespace internal
} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_LINE_ARC_2_H

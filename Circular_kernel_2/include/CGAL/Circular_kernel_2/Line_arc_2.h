// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_LINE_ARC_2_H
#define CGAL_CIRCULAR_KERNEL_LINE_ARC_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/global_functions_circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles/internal_functions_on_roots_and_polynomial_1_2_and_2_2.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_line_2.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_line_arc_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Circular_kernel_2/Circular_arc_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>


namespace CGAL {
namespace internal {

template <class CK >
class Line_arc_2_base
{
public:
  typedef typename CK::FT                        FT;
  typedef typename CK::RT                        RT;
  typedef typename CK::Point_2                   Point_2;
  typedef typename CK::Line_2                    Line_2;
  typedef typename CK::Circle_2                  Circle_2;
  typedef typename CK::Circular_arc_2            Circular_arc_2;
  typedef typename CK::Circular_arc_point_2      Circular_arc_point_2;
  typedef typename CK::Root_of_2                 Root_of_2;
  typedef typename CK::Segment_2                 Segment_2;

private:
  typedef struct bit_field {
    unsigned char begin_less_xy_than_end:2;
  } bit_field;

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

    typedef std::vector<typename CK2_Intersection_traits<CK, Line_2, Circle_2>::type>
      solutions_container;

    solutions_container solutions;
    CGAL::CircularFunctors::intersect_2<CK>
      ( l, c, std::back_inserter(solutions) );
    typename solutions_container::iterator it = solutions.begin();

    CGAL_kernel_precondition( it != solutions.end() );
    // the circles intersect

    const std::pair<typename CK::Circular_arc_point_2, unsigned>*
      result = CGAL::Intersections::internal::intersect_get<std::pair<typename CK::Circular_arc_point_2, unsigned> >(*it);
    // get must have succeeded
    if ( result->second == 2 ) // double solution
      return result->first;
    if (b) return result->first;
    ++it;
    result = CGAL::Intersections::internal::intersect_get<std::pair<typename CK::Circular_arc_point_2, unsigned> >(*it);
    return result->first;
  }


public:
  Line_arc_2_base()
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    : id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id())
#endif
  {}

  Line_arc_2_base(const Line_2 &support,
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


  Line_arc_2_base(const Line_2 &support,
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
    typename Intersection_traits<CK, Line_2, Line_2>::result_type
      v = CGAL::Intersections::internal::intersection(support, l1, CK());
    CGAL_assertion(bool(v));

    const Point_2 *pt = CGAL::Intersections::internal::intersect_get<Point_2>(v);
    CGAL_assertion(pt != nullptr);
    _begin = Circular_arc_point_2(*pt);
    v = CGAL::Intersections::internal::intersection(support, l2, CK());
    const Point_2 *pt2 = CGAL::Intersections::internal::intersect_get<Point_2>(v);
    CGAL_assertion(pt2 != nullptr);
    _end = Circular_arc_point_2(*pt2);
    reset_flags();
  }

  Line_arc_2_base(const Line_2 &support,
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

  Line_arc_2_base(const Segment_2 &s)
    :_support(s.supporting_line())
#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    ,id_of_my_supporting_line(Circular_arc_2::circle_table.get_new_id())
#endif
  {
    _begin = Circular_arc_point_2(s.source());
    _end = Circular_arc_point_2(s.target());
    reset_flags();
  }


  Line_arc_2_base(const Point_2 &p1,
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
                                            const Line_arc_2_base& l,
                                            T& res) {
    if(c.id_of_my_supporting_circle == 0) return false;
    return Circular_arc_2::circle_table.template find<T>(c.id_of_my_supporting_circle,
                                                         l.id_of_my_supporting_line,
                                                         res);
  }

  template < class T >
  static void put_intersection_circle_line(const Circular_arc_2& c,
                                           const Line_arc_2_base& l,
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

}; // end class Line_arc_2_base

template <class CB, typename Base_CK>
class Filtered_bbox_line_arc_2_base : public Base_CK::Line_arc_2 {
  typedef Filtered_bbox_line_arc_2_base<CB,Base_CK> Self;
  typedef typename Base_CK::Line_arc_2 P_arc;

public:

  typedef typename CB::Point_2 Point_2;
  typedef typename CB::Line_2 Line_2;
  typedef typename CB::Segment_2 Segment_2;
  typedef typename CB::Circle_2 Circle_2;
  typedef typename CB::Circular_arc_point_2 Circular_arc_point_2;

  Filtered_bbox_line_arc_2_base() : P_arc(), bb(nullptr) {}

  Filtered_bbox_line_arc_2_base(const P_arc& arc) : P_arc(arc), bb(nullptr) {}

  Filtered_bbox_line_arc_2_base(const Line_2 &support,
                                const Circle_2 &l1, const bool b_l1,
                                const Circle_2 &l2, const bool b_l2)
    : P_arc(support,l1,b_l1,l2,b_l2), bb(nullptr)
  {}


  Filtered_bbox_line_arc_2_base(const Line_2 &support,
                                const Line_2 &l1,
                                const Line_2 &l2)
    : P_arc(support,l1,l2), bb(nullptr)
  {}

  Filtered_bbox_line_arc_2_base(const Line_2 &support,
                                const Circular_arc_point_2 &begin,
                                const Circular_arc_point_2 &end)
    : P_arc(support, begin, end) , bb(nullptr)
  {}


  Filtered_bbox_line_arc_2_base(const Segment_2 &s)
    : P_arc(s) , bb(nullptr)
  {}


  Filtered_bbox_line_arc_2_base(const Point_2 &p1,
                                const Point_2 &p2)
    : P_arc(p1,p2) , bb(nullptr)
  {}


  Filtered_bbox_line_arc_2_base(const Filtered_bbox_line_arc_2_base &c)
    : P_arc(c), bb(c.bb ? new Bbox_2(*(c.bb)) : nullptr)
  {}

  Filtered_bbox_line_arc_2_base& operator=(const Self& c)
  {
    if(this != &c)
    {
      this->P_arc::operator=(c);

      if (bb != nullptr){
        delete bb;
      }
      bb = c.bb ? new Bbox_2(*(c.bb)) : nullptr;
    }
    return *this;
  }

  ~Filtered_bbox_line_arc_2_base() { if(bb) delete bb; }

  Bbox_2 bbox() const
  {
    if(bb==nullptr)
      bb=new Bbox_2(P_arc::bbox());
    return *bb;
  }

  bool has_no_bbox() const
  { return (bb==nullptr);}

private:

  mutable Bbox_2 *bb;

}; // end class Filtered_bbox_line_arc_2_base

/* template < typename CK > */
/*     std::ostream & */
/*     operator<<(std::ostream & os, const Line_arc_2_base<CK> &a) */
/*     { */

/*       return os << a.supporting_line() << " " */
/*                 << a.source() << " " */
/*                 << a.target() << " "; */
/*     } */

/*   template < typename CK > */
/*   std::istream & */
/*   operator>>(std::istream & is, Line_arc_2_base<CK> &a) */
/*   { */
/*     typename CK::Line_2 l; */
/*     typename CK::Circular_arc_point_2 p1; */
/*     typename CK::Circular_arc_point_2 p2; */
/*     is >> l >> p1 >> p2 ; */
/*     if (is) */
/*       a = Line_arc_2_base<CK>(l, p1, p2); */
/*     return is; */
/*   } */


} // namespace internal
} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_LINE_ARC_2_H

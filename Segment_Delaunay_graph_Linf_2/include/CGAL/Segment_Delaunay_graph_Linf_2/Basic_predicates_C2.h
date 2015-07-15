// Copyright (c) 2015  Università della Svizzera italiana.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_BASIC_PREDICATES_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_BASIC_PREDICATES_C2_H


#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_x_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_y_2.h>

#include <CGAL/Polychain_2.h>

#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

template<class K>
struct Basic_predicates_C2
 : public CGAL::SegmentDelaunayGraph_2::Basic_predicates_C2<K>
{
public:
  //-------------------------------------------------------------------
  // TYPES
  //-------------------------------------------------------------------

  typedef CGAL::SegmentDelaunayGraph_2::Basic_predicates_C2<K> Base;

  typedef typename K::RT                  RT;
  typedef typename K::FT                  FT;
  typedef typename K::Point_2             Point_2;
  typedef typename K::Segment_2           Segment_2;
  typedef typename K::Site_2              Site_2;
  typedef typename K::Oriented_side       Oriented_side;
  typedef typename K::Bounded_side        Bounded_side;
  typedef typename K::Comparison_result   Comparison_result;
  typedef typename K::Sign                Sign;
  typedef typename K::Orientation         Orientation;
  typedef typename K::Compute_scalar_product_2
                        Compute_scalar_product_2;
  typedef typename K::Boolean             Boolean;
  typedef typename K::Direction_2         Direction_2;
  typedef typename K::Vector_2            Vector_2;
  typedef typename K::Compare_x_2         Compare_x_2;
  typedef typename K::Compare_y_2         Compare_y_2;

  typedef typename CGAL::Polychainline_2<K> Polychainline_2;
  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>  Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K>
          Are_same_segments_2;

  typedef SegmentDelaunayGraph_2::Compare_x_2<K> Compare_x_2_Sites_Type;
  typedef SegmentDelaunayGraph_2::Compare_y_2<K> Compare_y_2_Sites_Type;

  typedef unsigned int Bearing;

private:
  typedef typename K::Intersections_tag ITag;
  typedef Bisector_Linf<K>              Bisector_Linf_Type;

  typedef typename Algebraic_structure_traits<RT>::Algebraic_category RT_Category;
  typedef typename Algebraic_structure_traits<FT>::Algebraic_category FT_Category;
public:

  typedef typename Base::Line_2               Line_2;
  typedef typename Base::Homogeneous_point_2  Homogeneous_point_2;

public:    //    compute_supporting_line(q.supporting_segment(), a1, b1, c1);
    //    compute_supporting_line(r.supporting_segment(), a2, b2, c2);

  //-------------------------------------------------------------------
  // BASIC CONSTRUCTIONS
  //-------------------------------------------------------------------

  using Base::compute_supporting_line;

  // compute horizontal line l that goes through p,
  // and leaves q on the oriented side s
  // s: has to be either +1 or -1 (not 0)
  // q: should not be on line l
  static
  Line_2
  compute_horizontal_side_line(
      const Point_2& p, const Point_2& q, Oriented_side s)
  {
    CGAL_precondition(s != ON_ORIENTED_BOUNDARY);

    RT b, c;

    b = RT(1);
    c = - p.y();

    // direction is (1, 0)

    Compare_y_2 cmpy;
    if (((cmpy(q, p) == LARGER) &&
         (s == ON_NEGATIVE_SIDE)   ) ||
        ((cmpy(q, p) == SMALLER) &&
         (s == ON_POSITIVE_SIDE)   )   ) {
      b = -b;
      c = -c;
    }
    return Line_2(RT(0), b, c);
  }

  // compute vertical line l that goes through p,
  // and leaves q on the oriented side s
  // s: has to be either +1 or -1 (not 0)
  // q: should not be on line l
  static
  Line_2
  compute_vertical_side_line(
      const Point_2& p, const Point_2& q, Oriented_side s)
  {
    CGAL_precondition(s != ON_ORIENTED_BOUNDARY);

    RT a, c;

    a = RT(1);
    c = - p.x();

    // direction is (0, -1)

    Compare_x_2 cmpx;
    if (((cmpx(q, p) == LARGER) &&
         (s == ON_NEGATIVE_SIDE)   ) ||
        ((cmpx(q, p) == SMALLER) &&
         (s == ON_POSITIVE_SIDE)   )   ) {
      a = -a;
      c = -c;
    }
    return Line_2(a, RT(0), c);
  }


  //using Base::compute_projection;

  static
  Homogeneous_point_2
  compute_linf_projection_hom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( ! l.is_degenerate() );
    CGAL_precondition(
        (CGAL::sign(l.a()) != ZERO) || (CGAL::sign(l.b()) != ZERO) );

    Sign signa = CGAL::sign(l.a());
    Sign signb = CGAL::sign(l.b());

    RT hx, hy, hw;

    if (signa == ZERO) { // l is horizontal
      // l.a() == 0  =>  y = -c/b
      hx = p.x() * l.b();
      hy = - l.c();
      hw = l.b();
    } else if (signb == ZERO) { // l is vertical
      // l.b() == 0  =>  x = -c/a
      hx = - l.c();
      hy = p.y() * l.a();
      hw = l.a();
    } else {
      // here both l.a() and l.b() are non-zero
      if ( signa == signb ) {
        hx = l.b() * ( p.x() - p.y() ) - l.c();
        hy = l.a() * ( p.y() - p.x() ) - l.c();
        hw = l.a() + l.b();
      } else { // signa != signb
        hx = -l.b() * ( p.x() + p.y() ) - l.c();
        hy = l.a() * ( p.x() + p.y() ) + l.c();
        hw = l.a() - l.b();
      }
    }

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Point_2
  compute_linf_projection_nonhom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( ! l.is_degenerate() );
    CGAL_precondition(
        (CGAL::sign(l.a()) != ZERO) || (CGAL::sign(l.b()) != ZERO) );

    Sign signa = CGAL::sign(l.a());
    Sign signb = CGAL::sign(l.b());

    RT hx, hy, hw;

    if (signa == ZERO) { // l is horizontal
      // l.a() == 0  =>  y = -c/b
      hx = p.x() * l.b();
      hy = - l.c();
      hw = l.b();
    } else if (signb == ZERO) { // l is vertical
      // l.b() == 0  =>  x = -c/a
      hx = - l.c();
      hy = p.y() * l.a();
      hw = l.a();
    } else {
      // here both l.a() and l.b() are non-zero
      if ( signa == signb ) {
        hx = l.b() * ( p.x() - p.y() ) - l.c();
        hy = l.a() * ( p.y() - p.x() ) - l.c();
        hw = l.a() + l.b();
      } else { // signa != signb
        hx = -l.b() * ( p.x() + p.y() ) - l.c();
        hy = l.a() * ( p.x() + p.y() ) + l.c();
        hw = l.a() - l.b();
      }
    }

    return Point_2(hx, hy, hw);
  }

  static
  Homogeneous_point_2
  compute_horizontal_projection_hom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( ! l.is_horizontal() );

    CGAL_precondition(
        (CGAL::sign(l.a()) != ZERO) );

    RT hx, hy, hw;

    hx = - l.b() * p.y() - l.c();
    hy = p.y() * l.a();
    hw = l.a();

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Point_2
  compute_horizontal_projection(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( ! l.is_horizontal() );
    CGAL_precondition(
        (CGAL::sign(l.a()) != ZERO) );

    RT hx, hy, hw;

    hx = - l.b() * p.y() - l.c();
    hy = p.y() * l.a();
    hw = l.a();

    return Point_2(hx, hy, hw);
  }

  static
  Homogeneous_point_2
  compute_vertical_projection_hom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( ! l.is_horizontal() );
    CGAL_precondition(
        (CGAL::sign(l.b()) != ZERO) );

    RT hx, hy, hw;

    hx = p.x() * l.b();
    hy = - l.a() * p.x() - l.c();
    hw = l.b();

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Point_2
  compute_vertical_projection(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( ! l.is_horizontal() );
    CGAL_precondition(
        (CGAL::sign(l.b()) != ZERO) );

    RT hx, hy, hw;

    hx = p.x() * l.b();
    hy = - l.a() * p.x() - l.c();
    hw = l.b();

    return Point_2(hx, hy, hw);
  }


  //using Base::projection_on_line;

  using Base::midpoint;

  using Base::compute_perpendicular;

  // compute_cw_perpendicular is opposite of compute_perpendicular
  static
  Line_2 compute_cw_perpendicular(const Line_2& l, const Point_2& p)
  {
    RT a, b, c;
    a = l.b();
    b = -l.a();
    c = -l.b() * p.x() + l.a() * p.y();
    return Line_2(a, b, c);
  }

  static
  Line_2 compute_linf_perpendicular(const Line_2& l, const Point_2& p)
  {
    RT a, b, c;
    a = RT( - CGAL::sign(l.b()) );
    b = RT( CGAL::sign(l.a()) );
    c = - a * p.x() - b * p.y();
    return Line_2(a, b, c);
  }

  using Base::opposite_line;

  // philaris: similar to compute_supporting_line
  static
  Line_2 compute_line_from_to(const Point_2& p, const Point_2&q)
  {
    RT a, b, c;
    a = p.y() - q.y();
    b = q.x() - p.x();

    CGAL_assertion((CGAL::sign(a) != ZERO) ||
                   (CGAL::sign(b) != ZERO))   ;

    c = p.x() * q.y() - q.x() * p.y();

    return Line_2(a, b, c);
  }

  // compute line from a point and a direction
  inline static
  Line_2 compute_line_dir(
      const Homogeneous_point_2& p, const Direction_2& d)
  {
    return Line_2( -d.dy()*p.hw(), d.dx()*p.hw(),
                   -(-d.dy()*p.hx() +d.dx()*p.hy()) );
  }

  // compute bisector of two parallel lines
  inline static
  Line_2 parallel_bis(const Line_2& lp, const Line_2& lq)
  {
    RT bisa, bisb, bisc;
    if ( CGAL::sign ( lq.a() ) != ZERO ) {
      bisa = RT(2) * lp.a() * lq.a();
      bisb = RT(2) * lp.a() * lq.b();
      bisc = lp.a() * lq.c() + lp.c() * lq.a();
    } else {
      bisa = RT(2) * lp.a() * lq.b();
      bisb = RT(2) * lp.b() * lq.b();
      bisc = lp.c() * lq.b() + lp.b() * lq.c();
    }
    return Line_2(bisa, bisb, bisc);
  }

  /* use point p for y coordinate of line */
  static
  Line_2 compute_horizontal_line_from_to(const Point_2& p, const Point_2&q)
  {
    RT b, c;
    Compare_x_2 cmpx;
    Comparison_result cmpxqp = cmpx(q,p);
    CGAL_assertion(cmpxqp != EQUAL);
    b = (cmpxqp == SMALLER) ? RT(-1) : RT(1);
    c = (cmpxqp == SMALLER) ? p.y() : -p.y();
    return Line_2(RT(0), b, c);
  }

  /* use point p for x coordinate of line */
  static
  Line_2 compute_vertical_line_from_to(const Point_2& p, const Point_2&q)
  {
    RT a, c;
    Compare_y_2 cmpy;
    Comparison_result cmpypq = cmpy(p,q);
    CGAL_assertion(cmpypq != EQUAL);
    a = (cmpypq == SMALLER) ? RT(-1) : RT(1);
    //a = RT(CGAL::sign(p.y() - q.y()));
    c = (cmpypq == SMALLER) ? p.x() : -p.x();
    return Line_2(a, RT(0), c);
  }


  // neg slope 45 degree line passing through p
  inline
  static
  Line_2 compute_neg_45_line_at(const Point_2 & p)
  {
    return Line_2(p.hw() , p.hw(), -p.hx()-p.hy());
  }

  // pos slope 45 degree line passing through p
  inline
  static
  Line_2 compute_pos_45_line_at(const Point_2 & p)
  {
    return Line_2(RT(1),RT(-1),p.y()-p.x());
  }

  // horizontal line passing through p
  inline
  static
  Line_2 compute_hor_line_at(const Point_2 & p)
  {
    return Line_2(RT(0), p.hw(), - p.hy());
  }

  // vertical line passing through p
  inline
  static
  Line_2 compute_ver_line_at(const Point_2 & p)
  {
    return Line_2(p.hw(), RT(0), - p.hx());
  }

  static
  RT compute_linf_distance(const Point_2& p, const Point_2& q)
  {
    return (CGAL::max)(
              CGAL::abs(p.x() - q.x()),
              CGAL::abs(p.y() - q.y()));
  }

  static
  std::pair<RT,RT>
  compute_linf_distance(const Point_2& p, const Line_2& l)
  {
    const RT nomin = CGAL::abs(l.a() * p.x() + l.b() * p.y() + l.c());
    const RT denom = CGAL::abs(
          l.a() +
          ( CGAL::sign(l.a()) == CGAL::sign(l.b())? l.b() : -l.b() ) );
    return std::pair<RT,RT>(nomin, denom);
  }

  static
  void compute_intersection_of_lines(
      const Line_2& l1, const Line_2& l2,
      RT& hx, RT& hy, RT& hw)
  {
    hx = l1.b() * l2.c() - l1.c() * l2.b();
    hy = l1.c() * l2.a() - l1.a() * l2.c();
    hw = l1.a() * l2.b() - l1.b() * l2.a();
  }

  inline
  static
  RT coord_at(const Line_2 &l, const RT & val, const bool return_y_coord)
  {
    return (return_y_coord)?
      (l.a() * val + l.c()) / (- l.b()) :
      (l.b() * val + l.c()) / (- l.a()) ;
  }

  inline
  static
  bool touch_same_side(
      const Site_2 & p, const Site_2 & q, const Line_2 & l,
      const bool samexpq, const bool pos_slope)
  {
    const RT common_coord = (samexpq) ? p.point().x() : p.point().y();
    const RT otherp = (samexpq) ? p.point().y() : p.point().x();
    const RT otherq = (samexpq) ? q.point().y() : q.point().x();
    const RT lineval = coord_at(l, common_coord, samexpq);
    return (CGAL::sign(lineval - otherp) == CGAL::sign(otherp - otherq)) ?
      samexpq == pos_slope :
      samexpq != pos_slope ;
  }

  inline
  static
  bool is_orth_dist_smaller_than_pt_dist(
      const RT & closest_coord, const Line_2 & l,
      const Site_2 & p, const Site_2 & q,
      const bool samexpq)
  {
    const RT lineval = coord_at(l, closest_coord, ! samexpq);
    return CGAL::abs(lineval - ((samexpq) ? p.point().x() : p.point().y()))
           <
           CGAL::abs((samexpq)? p.point().y() - q.point().y() :
                                p.point().x() - q.point().x()  ) ;
  }


public:
  //-------------------------------------------------------------------
  // BASIC PREDICATES
  //-------------------------------------------------------------------
  static
  Comparison_result
  compare_linf_distances_to_line(const Line_2& l, const Point_2& p,
                                    const Point_2& q)
  {
    Homogeneous_point_2 hp = compute_linf_projection_hom(l, p);
    Homogeneous_point_2 hq = compute_linf_projection_hom(l, q);

    RT dlp = (CGAL::max)(CGAL::abs(p.x() - hp.x()),
                         CGAL::abs(p.y() - hp.y()));

    RT dlq = (CGAL::max)(CGAL::abs(q.x() - hq.x()),
                         CGAL::abs(q.y() - hq.y()));

    Comparison_result crude = CGAL::compare(dlp, dlq);

    if (crude != EQUAL) {
      return crude;
    } else {
      CGAL_SDG_DEBUG(std::cout
          << "compare_linf_distances_to_line refining"
          << std::endl;);
      return crude;
    }
  }

  static
  Comparison_result
  compare_linf_distances_to_lines(const Point_2& p,
				     const Line_2& l1,
                                     const Line_2& l2)
  {
    Homogeneous_point_2 hl1 = compute_linf_projection_hom(l1, p);
    Homogeneous_point_2 hl2 = compute_linf_projection_hom(l2, p);

    RT dl1p = (CGAL::max)(CGAL::abs(hl1.x() - p.x()),
                          CGAL::abs(hl1.y() - p.y()));

    RT dl2p = (CGAL::max)(CGAL::abs(hl2.x() - p.x()),
                          CGAL::abs(hl2.y() - p.y()));

    Comparison_result crude = CGAL::compare(dl1p, dl2p);

    if (crude != EQUAL) {
      return crude;
    } else {
      CGAL_SDG_DEBUG(std::cout << "compare_linf_distances_to_lines refining"
                     << std::endl;);
      return crude;
    }
  }

  using Base::oriented_side_of_line;

  using Base::is_on_positive_halfspace;

  /* compares the Linf distances pq and pr;
     if the Linf distances pq are the same, then try to break ties
     by comparing the minima of the dx and dy distance components
   */
  static
  Comparison_result
  compare_distance_to_point_linf(
      const Point_2& p, const Point_2& q, const Point_2& r)
  {
    const RT pqdx = CGAL::abs(p.x()-q.x());
    const RT pqdy = CGAL::abs(p.y()-q.y());
    const bool pqdxlarger = CGAL::compare(pqdx, pqdy) == LARGER;
    const RT & pqmax = pqdxlarger ? pqdx : pqdy;
    const RT & pqmin = pqdxlarger ? pqdy : pqdx;

    const RT prdx = CGAL::abs(p.x()-r.x());
    const RT prdy = CGAL::abs(p.y()-r.y());
    const bool prdxlarger = CGAL::compare(prdx, prdy) == LARGER;
    const RT & prmax = prdxlarger ? prdx : prdy;
    const RT & prmin = prdxlarger ? prdy : prdx;

    const Comparison_result resmax = CGAL::compare(pqmax, prmax);

    if (resmax == EQUAL) {
      CGAL_SDG_DEBUG(std::cout <<
          "debug cmpdistlinf break ties with min" << std::endl;);
      return CGAL::compare(pqmin, prmin);
    } else {
      return resmax;
    }
  }

  static
  Bounded_side
  bounded_side_of_bbox(
      const Point_2& p, const Point_2& q, const Point_2& r)
  {
    // precondition: p, q, r should be monotone points.
    // Predicate bounded_side_of_bbox (P_bbox) returns:
    //  0 if p = q,
    //  0 if r = p or r = q (ON_BOUNDARY),
    // -1 if r is strictly outside the bounding box of p,q
    //    (ON_UNBOUNDED_SIDE),
    // +1 if r is strictly inside the bounding box of p,q
    //    (ON_BOUNDED_SIDE).
    // If p and q are on the same vertical or horizontal
    // line but are not the same point, then the bounding
    // box of p and q degenerates to the line segment pq.

    CGAL_SDG_DEBUG(std::cout << "debug bounded_side_of_bbox (p q r) = ("
                   << p << ") (" << q << ") (" << r << ")" << std::endl; );

    if ((CGAL::compare(p.x(), q.x()) == EQUAL) &&
        (CGAL::compare(p.y(), q.y()) == EQUAL)    ) {
      return ON_BOUNDARY;
    }

    Comparison_result cmpxpr, cmpxrq, cmpypr, cmpyrq;

    cmpxpr = CGAL::compare(p.x(), r.x());
    cmpxrq = CGAL::compare(r.x(), q.x());
    cmpypr = CGAL::compare(p.y(), r.y());
    cmpyrq = CGAL::compare(r.y(), q.y());

    Comparison_result comp =
      CGAL::compare(cmpxpr*cmpxrq + cmpypr*cmpyrq, 0);

    CGAL_SDG_DEBUG(std::cout << "debug bounded_side_of_bbox returns ";);

    switch(comp) {
      case SMALLER:
        CGAL_SDG_DEBUG(std::cout << "ON_UNBOUNDED_SIDE" << std::endl;);
        return ON_UNBOUNDED_SIDE;
      case EQUAL:
        CGAL_SDG_DEBUG(std::cout << "ON_BOUNDARY" << std::endl;);
        return ON_BOUNDARY;
      case LARGER:
        CGAL_SDG_DEBUG(std::cout << "ON_BOUNDED_SIDE" << std::endl;);
        return ON_BOUNDED_SIDE;
      default:
        CGAL_SDG_DEBUG(std::cout << "error: should never reach here";);
        CGAL_assertion( false );
        return ON_BOUNDARY;
    }
  }


  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the positive halfplane of oriented line l
  static
  Boolean
  intersects_segment_positive_halfplane(
      const Site_2 & s,
      const Line_2 & l)
  {
    Segment_2 seg = s.segment();

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    CGAL_SDG_DEBUG(std::cout
      << "debug: intersects_segment_positive_halfplane "
      << "s=" << s
      << " l=" << l.a() << " " << l.b() << " " << l.c()
                   << std::endl;);

    Oriented_side oslsrc = oriented_side_of_line(l, ssrc);
    Oriented_side osltrg = oriented_side_of_line(l, strg);

      CGAL_SDG_DEBUG(std::cout
          << "debug: intersects_segment_positive_halfplane "
          << "oslsrc=" << oslsrc << " osltrg=" << osltrg
                     << std::endl;);

    if ((oslsrc == ON_POSITIVE_SIDE) ||
        (osltrg == ON_POSITIVE_SIDE)   )
    {
      return true;
    } else {
      return false;
    }
  }


  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the negative halfplane of oriented line l
  static
  Boolean
  intersects_segment_negative_halfplane(
      const Site_2 & s,
      const Line_2 & l)
  {
    Segment_2 seg = s.segment();

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    CGAL_SDG_DEBUG(std::cout
        << "debug: intersects_segment_negative_halfplane "
        << "s=" << s
        << " l=" << l.a() << " " << l.b() << " " << l.c()
        << std::endl;);

    Oriented_side oslsrc = oriented_side_of_line(l, ssrc);
    Oriented_side osltrg = oriented_side_of_line(l, strg);

    CGAL_SDG_DEBUG(std::cout
        << "debug: intersects_segment_negative_halfplane "
        << "oslsrc=" << oslsrc << " osltrg=" << osltrg
        << std::endl;);

    if ((oslsrc == ON_NEGATIVE_SIDE) ||
        (osltrg == ON_NEGATIVE_SIDE)   )
    {
      return true;
    } else {
      return false;
    }
  }

  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // the only finite corner of the infinite box is corner
  // and if you traverse the infinite box ccw, then
  // you meet points in that order: q, corner, p
  static
  Boolean
  intersects_segment_interior_inf_box(const Site_2 & s,
      const Site_2 & q, const Site_2 & p,
      const Comparison_result & cmpxpq, const Comparison_result & cmpypq)
  {
    CGAL_assertion(cmpxpq != EQUAL);
    CGAL_assertion(cmpypq != EQUAL);
    CGAL_assertion(s.is_segment());
    const Segment_2 seg = s.segment();

    const Point_2 ssrc = seg.source();
    const Point_2 strg = seg.target();

    const Point_2 qq = q.point();
    const Point_2 pp = p.point();

    const bool eqcmp = cmpxpq == cmpypq;

    Are_same_points_2 same_points;
    Compare_x_2 cmpx;
    Compare_y_2 cmpy;

    bool is_ssrc_positive;
    if (same_points(q, s.source_site()) ||
        same_points(p, s.source_site())   ) {
      is_ssrc_positive = false;
    } else {
      const bool conflp = eqcmp ?
        (cmpx(pp, ssrc) == cmpxpq) : (cmpy(pp, ssrc) == cmpypq) ;
      const bool conflq = eqcmp ?
        (cmpy(ssrc, qq) == cmpypq) : (cmpx(ssrc, qq) == cmpxpq) ;
      is_ssrc_positive = (conflp && conflq);
    }
    if (is_ssrc_positive) {
      CGAL_SDG_DEBUG(std::cout << "debug is_segment_inside_inf_box "
                     << "src endpoint inside" << std::endl;);
      return true;
    }

    bool is_strg_positive;
    if (same_points(q, s.target_site()) ||
        same_points(p, s.target_site())   ) {
      is_strg_positive = false;
    } else {
      const bool conflp = eqcmp ?
        (cmpx(pp, strg) == cmpxpq) : (cmpy(pp, strg) == cmpypq) ;
      const bool conflq = eqcmp ?
        (cmpy(strg, qq) == cmpypq) : (cmpx(strg, qq) == cmpxpq) ;
      is_strg_positive = (conflp && conflq);
    }

    if (is_strg_positive) {
      CGAL_SDG_DEBUG(std::cout << "debug is_segment_inside_inf_box "
                     << "trg endpoint inside" << std::endl;);
      return true;
    } else {
      // here you have to check if the interior is inside

      CGAL_SDG_DEBUG(std::cout << "debug is_segment_inside_inf_box "
                     << "try for interior to be inside" << std::endl;);

      const Point_2 corner = eqcmp ?
        Point_2( pp.x(), qq.y() ) :
        Point_2( qq.x(), pp.y() ) ;

      // in fact, here you can intersect the segment
      // with the ray starting from corner and going to the
      // direction of the center of the infinite box

      const RT one(1);

      const Point_2 displaced ( corner.x() + (-cmpypq)*one ,
                                corner.y() + cmpxpq * one   );

      const Line_2 l = compute_line_from_to(corner, displaced);

      const Line_2 lseg = compute_supporting_line(s.supporting_site());

      RT hx, hy, hw;
      compute_intersection_of_lines(l, lseg, hx, hy, hw);

      if (CGAL::sign(hw) == ZERO) {
        return false;
      } else {
        const Point_2 ip ( hx, hy, hw );
        const Line_2 lqc = compute_line_from_to(qq, corner);
        const Line_2 lcp = compute_line_from_to(corner, pp);
        const Oriented_side os_lqc_ip = oriented_side_of_line(lqc, ip);
        const Oriented_side os_lcp_ip = oriented_side_of_line(lcp, ip);

        const Comparison_result cmpxsrcip = cmpx(ssrc, ip);
        const Comparison_result cmpysrcip = cmpy(ssrc, ip);
        const Comparison_result cmpxiptrg = cmpx(ip, strg);
        const Comparison_result cmpyiptrg = cmpy(ip, strg);

        // philaris: to check
        Boolean is_ip_inside_segment =
          (CGAL::sign(cmpxsrcip * cmpxiptrg +
                      cmpysrcip * cmpyiptrg   )) == POSITIVE;

        if ((os_lqc_ip == ON_POSITIVE_SIDE) &&
            (os_lcp_ip == ON_POSITIVE_SIDE) &&
            is_ip_inside_segment ) {
          return true;
        } else {
          return false;
        }
      }
    }
  } // end of intersects_segment_interior_inf_box

  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // the only finite corner of the infinite box is corner
  // and if you traverse the infinite box ccw, then
  // you meet points in that order: q, corner, p
  static
  Boolean
  intersects_segment_interior_inf_box(const Site_2 & s,
      const Site_2 & q, const Point_2 & corner,
      const Site_2 & p)
  {
    CGAL_assertion(s.is_segment());
    Segment_2 seg = s.segment();

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    Point_2 qq = q.point();
    Point_2 pp = p.point();

    Line_2 lqc = compute_line_from_to(qq, corner);
    Line_2 lcp = compute_line_from_to(corner, pp);

    Are_same_points_2 same_points;

    bool is_ssrc_positive;
    if (same_points(q, s.source_site()) ||
        same_points(p, s.source_site())   ) {
      is_ssrc_positive = false;
    } else {
      Oriented_side os_lqc_ssrc = oriented_side_of_line(lqc, ssrc);
      Oriented_side os_lcp_ssrc = oriented_side_of_line(lcp, ssrc);
      is_ssrc_positive =
        ((os_lqc_ssrc == ON_POSITIVE_SIDE) &&
         (os_lcp_ssrc == ON_POSITIVE_SIDE)    ) ;
    }

    bool is_strg_positive;
    if (same_points(q, s.target_site()) ||
        same_points(p, s.target_site())   ) {
      is_strg_positive = false;
    } else {
      Oriented_side os_lqc_strg = oriented_side_of_line(lqc, strg);
      Oriented_side os_lcp_strg = oriented_side_of_line(lcp, strg);
      is_strg_positive =
        ((os_lqc_strg == ON_POSITIVE_SIDE) &&
         (os_lcp_strg == ON_POSITIVE_SIDE)    ) ;
    }

    CGAL_SDG_DEBUG(std::cout << "debug qcp= (" << q << ") (" << corner
        << ") (" << p << ")"
        << " isssrcpos=" << is_ssrc_positive
        << " isstrgpos=" << is_strg_positive
        << std::endl;);

    if (is_ssrc_positive || is_strg_positive) {
      CGAL_SDG_DEBUG(std::cout << "debug is_segment_inside_inf_box "
                     << "endpoint inside" << std::endl;);
      return true;
    } else {
      // here you have to check if the interior is inside

      CGAL_SDG_DEBUG(std::cout << "debug is_segment_inside_inf_box "
                     << "try for interior to be inside" << std::endl;);

      // in fact, here you can intersect the segment
      // with the ray starting from corner and going to the
      // direction of the center of the infinite box

      Compare_x_2 cmpx;
      Compare_y_2 cmpy;

      Comparison_result cmpxpq = cmpx(pp,qq);
      Comparison_result cmpypq = cmpy(pp,qq);

      RT one(1);

      Point_2 displaced ( corner.x() + (-cmpypq)*one ,
                          corner.y() + cmpxpq * one   );

      Line_2 l = compute_line_from_to(corner, displaced);

      Line_2 lseg = compute_supporting_line(s.supporting_site());

      RT hx, hy, hw;

      compute_intersection_of_lines(l, lseg, hx, hy, hw);

      if (CGAL::sign(hw) == ZERO) {
        return false;
      } else {
        Point_2 ip ( hx/hw, hy/hw);
        Oriented_side os_lqc_ip = oriented_side_of_line(lqc, ip);
        Oriented_side os_lcp_ip = oriented_side_of_line(lcp, ip);

        Compare_x_2 cmpx;
        Compare_y_2 cmpy;

        Comparison_result cmpxsrcip = cmpx(ssrc, ip);
        Comparison_result cmpysrcip = cmpy(ssrc, ip);
        Comparison_result cmpxiptrg = cmpx(ip, strg);
        Comparison_result cmpyiptrg = cmpy(ip, strg);

        // philaris: to check
        Boolean is_ip_inside_segment =
          (CGAL::sign(cmpxsrcip * cmpxiptrg +
                      cmpysrcip * cmpyiptrg   )) == POSITIVE;

        if ((os_lqc_ip == ON_POSITIVE_SIDE) &&
            (os_lcp_ip == ON_POSITIVE_SIDE) &&
            is_ip_inside_segment ) {
          return true;
        } else {
          return false;
        }
      }
    }
  } // end of intersects_segment_interior_inf_box


  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // this infinite box is a 90 degree wedge defined
  // by the intersection of the halfplanes
  // with supporting lines lhor and lver, where
  // the halfplanes are both on the positive or negative
  // sides of the supporting lines
  static
  Boolean
  intersects_segment_side_of_wedge(const Site_2 & s,
      const Line_2 & lhor, const Line_2 & lver,
      Oriented_side orside)
  {
    CGAL_assertion(s.is_segment());
    Segment_2 seg = s.segment();

    CGAL_SDG_DEBUG(std::cout << "debug sofw s=" << s
                   << " orside=" << orside << std::endl;);

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    Oriented_side os_lhor_ssrc = oriented_side_of_line(lhor, ssrc);
    Oriented_side os_lver_ssrc = oriented_side_of_line(lver, ssrc);

    Oriented_side os_lhor_strg = oriented_side_of_line(lhor, strg);
    Oriented_side os_lver_strg = oriented_side_of_line(lver, strg);

    if (((os_lhor_ssrc == orside) &&
         (os_lver_ssrc == orside)) ||
        ((os_lhor_strg == orside) &&
         (os_lver_strg == orside))   ) {
          CGAL_SDG_DEBUG(std::cout
              << "debug intersects_segment_side_of_wedge "
              << "endpoint inside" << std::endl;);
      return true;
    } else {
      // here we have to check if the interior is inside

      CGAL_SDG_DEBUG(std::cout
          << "debug intersects_segment_side_of_wedge "
          << "try for interior to be inside" << std::endl;);

      // in fact, here you can intersect the segment
      // with the ray starting from corner and going to the
      // direction of the center of the infinite box

      // corner has homogenuous coordinates cx, cy, cw
      RT cx, cy, cw;
      compute_intersection_of_lines(lhor, lver, cx, cy, cw);

      CGAL_assertion( CGAL::sign(cw) != ZERO );

      Point_2 corner ( cx, cy, cw );

      CGAL_SDG_DEBUG(std::cout << "debug corner=" << corner << std::endl;);

      RT one(1);

      Point_2 displaced (
          corner.x() + ( (+orside)*CGAL::sign(lver.a()) ) * one ,
          corner.y() +   (+orside)*CGAL::sign(lhor.b())   * one   );

      CGAL_SDG_DEBUG(std::cout
          << "debug displaced=" << displaced << std::endl;);

      Line_2 l = compute_line_from_to(corner, displaced);

      Line_2 lseg = compute_supporting_line(s.supporting_site());

      RT hx, hy, hw;

      CGAL_SDG_DEBUG(std::cout
          << "debug: intersects_segment_side_of_wedge "
          << " l=" << l.a() << " " << l.b() << " " << l.c()
          << " lseg=" << lseg.a() << " " << lseg.b() << " " << lseg.c()
          << std::endl;);

      compute_intersection_of_lines(l, lseg, hx, hy, hw);

      if (CGAL::sign(hw) == ZERO) {
        CGAL_SDG_DEBUG(std::cout
            << "debug l and lseg are parallel" << std::endl;);
        return false;
      } else {
        Point_2 ip ( hx, hy, hw );
        CGAL_SDG_DEBUG(std::cout << "debug ip=" << ip << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug ip_hom="
            << hx << ' ' << hy << ' ' << hw << std::endl;);
        Oriented_side os_lhor_ip = oriented_side_of_line(lhor, ip);
        Oriented_side os_lver_ip = oriented_side_of_line(lver, ip);

        Compare_x_2 cmpx;
        Compare_y_2 cmpy;

        Comparison_result cmpxsrcip = cmpx(ssrc, ip);
        Comparison_result cmpysrcip = cmpy(ssrc, ip);
        Comparison_result cmpxiptrg = cmpx(ip, strg);
        Comparison_result cmpyiptrg = cmpy(ip, strg);

        // philaris: to check
        Boolean is_ip_inside_segment =
          (CGAL::sign(cmpxsrcip * cmpxiptrg +
                      cmpysrcip * cmpyiptrg   )) == POSITIVE;

        if ((os_lhor_ip == orside) &&
            (os_lver_ip == orside) &&
            is_ip_inside_segment      ) {
          return true;
        } else {
          return false;
        }
      }
    }
  } // end of intersects_segment_side_of_wedge

  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // this infinite box is a 90 degree wedge defined
  // by the intersection of the positive halfplanes
  // with supporting lines lhor and lver
  static
  Boolean
  intersects_segment_positive_of_wedge(const Site_2 & s,
      const Line_2 & lhor, const Line_2 & lver)
  {
    return intersects_segment_side_of_wedge(
        s, lhor, lver, ON_POSITIVE_SIDE);
  }

  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // this infinite box is a 90 degree wedge defined
  // by the intersection of the positive halfplanes
  // with supporting lines lhor and lver
  static
  Boolean
  intersects_segment_negative_of_wedge(const Site_2 & s,
      const Line_2 & lhor, const Line_2 & lver)
  {
    return intersects_segment_side_of_wedge(
        s, lhor, lver, ON_NEGATIVE_SIDE);
  }

  // returns true if and only if
  // the interior of t has non-empty intersection
  // with the interior of the following infinite box:
  // the only finite corner of the infinite box is p
  // and this infinite box is a 90 degree wedge defined
  // by the intersection of the halfplanes
  // with supporting lines lhor and lver through p, where
  // the halfplanes are both on the positive or negative
  // sides of the supporting lines. The positive or negative
  // side depends on the relative position of p with respect to s

  static
  Boolean
  intersects_segment_interior_inf_wedge_sp(const Site_2 & s,
                                           const Site_2 & p,
                                           const Site_2 & t)
  {
    CGAL_assertion( t.is_segment() );
    CGAL_assertion( s.is_segment() );
    CGAL_assertion(! is_site_h_or_v(s));

    Segment_2 seg = s.segment();

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    Point_2 pp = p.point();

    Sign dxs = CGAL::sign(strg.x() - ssrc.x());
    Sign dys = CGAL::sign(strg.y() - ssrc.y());

    Line_2 lseg = compute_supporting_line(s.supporting_site());
    Oriented_side os_lseg_p = oriented_side_of_line(lseg, pp);

    CGAL_assertion( os_lseg_p != ON_ORIENTED_BOUNDARY );

    //sandeep: lhor, lver remains same for left turn and right turn

    if (dxs == NEGATIVE && dys == NEGATIVE) {

      Line_2 lhor = Line_2(0,1,-pp.y());
      Line_2 lver = Line_2(-1,0,pp.x());

      return intersects_segment_side_of_wedge(t,
                                       lhor, lver,
                                       os_lseg_p);
    } else if (dxs == POSITIVE && dys == NEGATIVE) {

      Line_2 lhor = Line_2(0,-1,pp.y());
      Line_2 lver = Line_2(-1,0,pp.x());

      return intersects_segment_side_of_wedge(t,
                                       lhor, lver,
                                       os_lseg_p);
    } else if (dxs == POSITIVE && dys == POSITIVE) {

      Line_2 lhor = Line_2(0,-1,pp.y());
      Line_2 lver = Line_2(1,0,-pp.x());

      return intersects_segment_side_of_wedge(t,
                                       lhor, lver,
                                       os_lseg_p);
    } else {//dxs == NEGATIVE and dys == POSITIVE

      Line_2 lhor = Line_2(0,1,-pp.y());
      Line_2 lver = Line_2(1,0,-pp.x());

      return intersects_segment_side_of_wedge(t,
                                       lhor, lver,
                                       os_lseg_p);
    }

  } // end of intersects_segment_interior_inf_wedge_sp


  // returns true if and only if
  // the interior of s has non-empty intersection
  // with the interior of the bounding box of q, p
  // precondition: the bounding box should be non-trivial,
  // i.e., it should not be a segment
  static
  Boolean
  intersects_segment_interior_bbox(const Site_2 & s,
      const Site_2 & q,
      const Site_2 & p)
  {
    CGAL_precondition(s.is_segment());
    CGAL_precondition(p.is_point());
    CGAL_precondition(q.is_point());

    Point_2 pp = p.point();
    Point_2 qq = q.point();

    CGAL_assertion_code( Compare_x_2 cmpx; )
    CGAL_assertion_code( Compare_y_2 cmpy; )
    CGAL_assertion(cmpx(pp,qq) != EQUAL);
    CGAL_assertion(cmpy(pp,qq) != EQUAL);

    Point_2 corner1 ( pp.x(), qq.y());
    Point_2 corner2 ( qq.x(), pp.y());

    if (CGAL::orientation( qq, corner1, pp ) == LEFT_TURN) {
      return intersects_segment_interior_inf_box(s, q, corner1, p)
          && intersects_segment_interior_inf_box(s, p, corner2, q);
    } else {
      return intersects_segment_interior_inf_box(s, q, corner2, p)
          && intersects_segment_interior_inf_box(s, p, corner1, q);
    }
  } // end of intersects_segment_interior_bbox

  // returns true if and only if
  // the non-horizontal/non-vertical segment at site s
  // has positive slope
  static
  Boolean
  has_positive_slope(const Site_2 & s)
  {
    CGAL_precondition(s.is_segment());
    CGAL_precondition(! is_site_h_or_v(s));
    Compare_x_2 cmpx;
    Compare_y_2 cmpy;
    Point_2 src = s.supporting_site().source();
    Point_2 trg = s.supporting_site().target();
    return cmpx(src, trg) == cmpy(src, trg);
  }

  inline
  static
  bool
  has_positive_slope(const Line_2 & l) {
    return (CGAL::sign(l.a()) + CGAL::sign(l.b()) == 0);
  }

  inline
  static
  Boolean
  have_same_slope(const Site_2 & s, const Site_2 & t)
  {
    CGAL_precondition(s.is_segment());
    CGAL_precondition(t.is_segment());
    Compare_x_2 cmpx;
    Compare_y_2 cmpy;
    Point_2 ssrc = s.supporting_site().source();
    Point_2 strg = s.supporting_site().target();
    Comparison_result scmpx = cmpx(ssrc, strg);
    Comparison_result scmpy = cmpy(ssrc, strg);
    Point_2 tsrc = t.supporting_site().source();
    Point_2 ttrg = t.supporting_site().target();
    Comparison_result tcmpx = cmpx(tsrc, ttrg);
    Comparison_result tcmpy = cmpy(tsrc, ttrg);
    CGAL_SDG_DEBUG(std::cout << "debug have_same_slope"
        << " scmpx=" << scmpx << " scmpy=" << scmpy
        << " tcmpx=" << tcmpx << " tcmpy=" << tcmpy
        << std::endl;);
    if (   ((scmpx == EQUAL) && (tcmpx == EQUAL)) // vertical
        || ((scmpy == EQUAL) && (tcmpy == EQUAL)) // horizontal
        || ((scmpx == scmpy) && (tcmpx == tcmpy)) // positive
        || ((scmpx != EQUAL) && (scmpy != EQUAL) &&
            (tcmpx != EQUAL) && (tcmpy != EQUAL) &&
            (scmpx != scmpy) && (tcmpx != tcmpy)) // negative
       ) {
      return true;
    } else {
      return false;
    }
  }

  inline
  static
  Boolean
  is_site_h_or_v(const Site_2 & s)
  {
    return is_site_horizontal(s) || is_site_vertical(s);
  }

  inline
  static
  Boolean
  is_site_horizontal(const Site_2 & s)
  {
    CGAL_assertion(s.is_segment());
    return s.supporting_site().segment().is_horizontal();
  }

  inline
  static
  Boolean
  is_site_vertical(const Site_2 & s)
  {
    CGAL_assertion(s.is_segment());
    return s.supporting_site().segment().is_vertical();
  }

  inline
  static
  Boolean
  is_line_h_or_v(const Line_2 & l)
  {
    return (CGAL::sign(l.a()) == ZERO) || (CGAL::sign(l.b()) == ZERO);
  }

  inline
  static
  Boolean
  test_star(const Site_2 & p, const Site_2 & u,
            const Site_2 & v, const Site_2 & t) {
    CGAL_precondition(p.is_point());
    CGAL_precondition(u.is_segment());
    CGAL_precondition(v.is_segment());
    CGAL_precondition(t.is_segment());

    Are_same_points_2 same_points;

    Point_2 pu =
      (same_points(p, u.source_site()) ? u.target_site() : u.source_site())
       .point();
    Point_2 pv =
      (same_points(p, v.source_site()) ? v.target_site() : v.source_site())
       .point();
    Point_2 pt =
      (same_points(p, t.source_site()) ? t.target_site() : t.source_site())
       .point();

    Orientation oupt = CGAL::orientation(pu, p.point(), pt);
    Orientation otpv = CGAL::orientation(pt, p.point(), pv);

    return (oupt == LEFT_TURN) && (otpv == LEFT_TURN);
  }

  static
  Boolean
  are_in_same_open_halfspace_of(
      const Site_2 & p, const Site_2 & q, const Site_2 & r)
  {
    CGAL_precondition(p.is_point() && q.is_point() && r.is_segment());
    Line_2 lseg = compute_supporting_line(r.supporting_site());
    Oriented_side os_lseg_p = oriented_side_of_line(lseg, p.point());
    if ( os_lseg_p == ON_ORIENTED_BOUNDARY ) {
      return false;
    }
    Oriented_side os_lseg_q = oriented_side_of_line(lseg, q.point());
    return os_lseg_p == os_lseg_q;
  }

  inline
  static
  RT horseg_y_coord(const Site_2 & s) {
    CGAL_assertion(s.is_segment());
    CGAL_assertion(is_site_horizontal(s));
    return s.supporting_site().source_site().point().y();
  }

  inline
  static
  RT verseg_x_coord(const Site_2 & s) {
    CGAL_assertion(s.is_segment());
    CGAL_assertion(is_site_vertical(s));
    return s.supporting_site().source_site().point().x();
  }

  inline
  static
  RT hvseg_coord(const Site_2 & s, const bool is_hor) {
    CGAL_assertion(s.is_segment());
    CGAL_assertion(is_site_horizontal(s) == is_hor);
    CGAL_assertion(is_site_vertical(s) == (! is_hor));
    return is_hor ? s.supporting_site().source_site().point().y() :
                    s.supporting_site().source_site().point().x() ;
  }

  inline static
  bool check_if_exact(const Site_2& , unsigned int ,
		      const Tag_false&)
  {
    return true;
  }

  inline static
  bool check_if_exact(const Site_2& s, unsigned int i,
		      const Tag_true&)
  {
    return s.is_input(i);
  }

  // determines if the segment s is on the positive halfspace as
  // defined by the supporting line of the segment supp; the line l
  // is supposed to be the supporting line of the segment supp and we
  // pass it so that we do not have to recompute it
  static bool
  is_on_positive_halfspace(const Site_2& supp,
			   const Site_2& s, const Line_2& l)
  {
    CGAL_precondition( supp.is_segment() && s.is_segment() );
    Are_same_points_2 same_points;
    Are_same_segments_2 same_segments;

    if ( same_segments(supp.supporting_site(),
		       s.supporting_site()) ) {
      return false;
    }

    if ( same_points(supp.source_site(), s.source_site()) ||
	 same_points(supp.target_site(), s.source_site()) ) {
      return oriented_side_of_line(l, s.target()) == ON_POSITIVE_SIDE;
    }

    if ( same_points(supp.source_site(), s.target_site()) ||
	 same_points(supp.target_site(), s.target_site()) ) {
      return oriented_side_of_line(l, s.source()) == ON_POSITIVE_SIDE;
    }

    ITag itag;

    if ( !check_if_exact(s, 0, itag) &&
	 same_segments(supp.supporting_site(),
		       s.crossing_site(0)) ) {
      return oriented_side_of_line(l, s.target()) == ON_POSITIVE_SIDE;
    }

    if ( !check_if_exact(s, 1, itag) &&
	 same_segments(supp.supporting_site(),
		       s.crossing_site(1)) ) {
      return oriented_side_of_line(l, s.source()) == ON_POSITIVE_SIDE;
    }

    return Base::is_on_positive_halfspace(l, s.segment());
  }

  // The bearing of a line l is defined as a number in {0, 1, ..., 7},
  // encoding each permissible ordered pair (sign(l.a), sign(l.b)).
  // There are 8 permissible pairs; pair (ZERO, ZERO) is not allowed.
  // A line with pair (NEG, POS) has bearing 0 and as this line rotates
  // counterclockwise we get consecutive bearings. Observe that
  // axis-parallel lines have odd bearings (1, 3, 5, 7).

  inline static Bearing
  bearing(const Line_2 & l) {
    const Sign sa = CGAL::sign(l.a());
    const Sign sb = CGAL::sign(l.b());
    if (sa == NEGATIVE) {
      return 1-sb;
    } else if (sa == ZERO) {
      return (sb == NEGATIVE) ? 3 : 7;
    } else { // sa == POSITIVE
      return 5+sb;
    }
  }

private:
  inline static bool
  have_same_bearing(const Line_2 & l1, const Line_2 & l2) {
    return (CGAL::sign(l1.a()) == CGAL::sign(l2.a())) &&
           (CGAL::sign(l1.b()) == CGAL::sign(l2.b()))    ;
  }

  inline static bool
  have_opposite_bearing(const Line_2 & l1, const Line_2 & l2) {
    return (CGAL::sign(l1.a()) == -CGAL::sign(l2.a())) &&
           (CGAL::sign(l1.b()) == -CGAL::sign(l2.b()))    ;
  }

  inline static bool
  bearing_outside(
      const Bearing bprev,
      const Bearing br,
      const Bearing bnext)
  {
    CGAL_assertion(bprev <= 7);
    CGAL_assertion(br <= 7);
    CGAL_assertion(bnext <= 7);
    CGAL_assertion(bprev != bnext);
    CGAL_assertion(bprev != br);
    CGAL_assertion(bnext != br);
    Bearing i(bprev);
    while (true) {
      i = (i + 1) % 8;
      if (i == bnext) {
        return true;
      } else if (i == br) {
        return false;
      }
    }
  }

public:
  // absolute bearing difference between low and high
  inline static unsigned int
  bearing_diff(const Bearing low, const Bearing high)
  {
    CGAL_assertion(low <= 7);
    CGAL_assertion(high <= 7);
    return high > low ? high - low : 8 + high - low;
  }

  // Orient the segments p, q, r so that they go counterclockwise
  // around the Linf square. The orientations are saved as lines
  // l[0], l[1], l[2] in the l array.
  static void
  orient_lines_linf(const Site_2& p, const Site_2& q, const Site_2& r,
      Line_2 l[])
  {
    CGAL_precondition( p.is_segment() && q.is_segment() &&
		       r.is_segment() );

    l[0] = compute_supporting_line(p.supporting_site());
    l[1] = compute_supporting_line(q.supporting_site());
    l[2] = compute_supporting_line(r.supporting_site());

    bool is_oriented[3] = {false, false, false};

    if ( is_on_positive_halfspace(p, q, l[0]) ||
	 is_on_positive_halfspace(p, r, l[0]) ) {
      is_oriented[0] = true;
    } else {
      l[0] = opposite_line(l[0]);
      if ( is_on_positive_halfspace(p, q, l[0]) ||
	   is_on_positive_halfspace(p, r, l[0]) ) {
	is_oriented[0] = true;
      } else {
	l[0] = opposite_line(l[0]);
      }
    }

    if ( is_on_positive_halfspace(q, p, l[1]) ||
	 is_on_positive_halfspace(q, r, l[1]) ) {
      is_oriented[1] = true;
    } else {
       l[1] = opposite_line(l[1]);
      if ( is_on_positive_halfspace(q, p, l[1]) ||
	   is_on_positive_halfspace(q, r, l[1]) ) {
	is_oriented[1] = true;
      } else {
	l[1] = opposite_line(l[1]);
      }
    }

    if ( is_on_positive_halfspace(r, p, l[2]) ||
	 is_on_positive_halfspace(r, q, l[2]) ) {
      is_oriented[2] = true;
    } else {
      l[2] = opposite_line(l[2]);
      if ( is_on_positive_halfspace(r, p, l[2]) ||
	   is_on_positive_halfspace(r, q, l[2]) ) {
	is_oriented[2] = true;
      } else {
	l[2] = opposite_line(l[2]);
      }
    }

    if ( is_oriented[0] && is_oriented[1] && is_oriented[2] ) {
      return;
    }

    int i_no(-1);
    for (int i = 0; i < 3; i++) {
      if ( !is_oriented[i] ) {
	i_no = i;
	CGAL_assertion( is_oriented[(i+1)%3] && is_oriented[(i+2)%3] );
	break;
      }
    }

    CGAL_assertion( i_no != -1 );

    CGAL_SDG_DEBUG( std::cout << "orient_lines_linf i_no: "
        << p << ' ' << q << ' ' << ' ' << r << " i_no=" << i_no
        << std::endl;);

    for (int j = 1; j < 3; j++) {
      if (have_same_bearing(l[i_no], l[(i_no+j)%3])) {
        l[i_no] = opposite_line(l[i_no]);
        is_oriented[i_no] = true;
      } else if (have_opposite_bearing(l[i_no], l[(i_no+j)%3])) {
        is_oriented[i_no] = true;
      }
      if (is_oriented[i_no]) break;
    }

    if (is_oriented[i_no]) {
      return;
    }

    CGAL_SDG_DEBUG( std::cout << "orient_lines_linf lonely bearing: "
        << p << ' ' << q << ' ' << ' ' << r << " i_no=" << i_no
        << std::endl;);
    const unsigned int iprev = (i_no+2)%3;
    const unsigned int inext = (i_no+1)%3;
    const Bearing bprev = bearing(l[iprev]);
    const Bearing bnext = bearing(l[inext]);
    CGAL_SDG_DEBUG( std::cout << "orient_lines_linf"
        << " bprev=" << bprev << " bnext=" << bnext
        << std::endl;);
    CGAL_assertion(bprev != bnext);
    CGAL_assertion_code( const Bearing diffbear = (bnext - bprev)%8 );
    CGAL_assertion(diffbear != 1);
    CGAL_assertion(diffbear != 7);
    const Bearing br = bearing(l[i_no]);
    if ( bearing_outside(bprev, br, bnext) ) {
      CGAL_assertion_code( const Bearing bropp = (br+4)%8 );
      CGAL_assertion( ! bearing_outside(bprev, bropp, bnext) );
      l[i_no] = opposite_line(l[i_no]);
      is_oriented[i_no] = true;
    } else {
      CGAL_assertion( ! bearing_outside(bprev, br, bnext) );
      const Bearing bropp = (br+4)%8;
      if ( ! bearing_outside(bprev, bropp, bnext) ) {
        CGAL_assertion( diffbear == 6 );
        CGAL_assertion( br % 2 == 1 ); // undecided is axis-parallel
        const Site_2 & sprev = iprev == 0 ? p : (iprev == 1? q : r);
        Bearing brcorrect = (bprev+5)%8;
        if (Base::is_on_positive_halfspace(l[inext], sprev.segment())) {
          brcorrect = (bprev+1)%8;
        } else {
          CGAL_assertion_code(
              const Site_2 & snext = inext == 0 ? p : (inext == 1? q : r) );
          CGAL_assertion(
              Base::is_on_positive_halfspace(l[iprev], snext.segment()) );
        }
        if (brcorrect == bropp) {
          l[i_no] = opposite_line(l[i_no]);
        } else {
          CGAL_assertion(brcorrect == br);
        }
      }
      is_oriented[i_no] = true;
    }
    CGAL_assertion(is_oriented[i_no]);
  }

  inline static bool
  are_parallel_lines(const Line_2 & lp, const Line_2 & lq) {
    return lp.a() * lq.b() == lq.a() * lp.b();
  }

  inline static Direction_2
  direction(const Line_2 & l) {
    return Direction_2(l.b(), -l.a());
  }

  // compute correct bisector direction of two consecutive lines,
  // that are already correctly oriented from sites with orient_lines_linf
  static Direction_2
  dir_from_lines(const Line_2& lp, const Line_2& lq)
  {
    Bisector_Linf_Type linf_bisect_direction;
    const unsigned int bdiff = bearing_diff(bearing(lp), bearing(lq));
    CGAL_assertion(bdiff != 0);
    if (bdiff < 4) {
      return linf_bisect_direction(direction(lq), -direction(lp));
    } else if (bdiff > 4) {
      return linf_bisect_direction(direction(lp), -direction(lq));
    } else {
      CGAL_assertion(bdiff == 4);
      const Sign sgn = CGAL::sign(lp.a() * lq.b() - lq.a() * lp.b());
      if (sgn == POSITIVE) {
        return linf_bisect_direction(direction(lq), -direction(lp));
      } else {
        CGAL_assertion(sgn == NEGATIVE);
        return linf_bisect_direction(direction(lp), -direction(lq));
      }
    }
  }

  // Compute the bisecting line between sites p and q, such that
  // the sites are correctly oriented according to lines lp and lq,
  // respectively
  static Line_2
  bisector_linf_line(const Site_2& p, const Site_2& q,
      const Line_2 & lp, const Line_2 & lq)
  {
    if (are_parallel_lines(lp, lq)) {
      return parallel_bis(lp, lq);
    } else {
      return bisector_linf_line_nonpar(p, q, lp, lq);
    }
  }

  inline static Line_2
  bisector_linf_line_nonpar(const Site_2& p, const Site_2& q,
      const Line_2 & lp, const Line_2 & lq)
  {
    const bool is_psrc_q = is_endpoint_of(p.source_site(), q);
    const bool is_ptrg_q = is_endpoint_of(p.target_site(), q);
    const bool have_common_pq = is_psrc_q || is_ptrg_q;
    Homogeneous_point_2 xpq;
    if (have_common_pq) {
      xpq = is_psrc_q ? p.source() : p.target();
    } else {
      RT hx, hy, hw;
      compute_intersection_of_lines(lp, lq, hx, hy, hw);
      CGAL_SDG_DEBUG( std::cout << "debug xpq hom="
        << hx << ' ' << hy << ' ' << hw << std::endl; );
      xpq = Homogeneous_point_2(hx, hy, hw);
    }
    const Direction_2 dirbpq = dir_from_lines(lp, lq);
    CGAL_SDG_DEBUG( std::cout << "debug xpq hom="
        << xpq.hx() << ' ' << xpq.hy() << ' ' << xpq.hw()
        << " dirbpq=" << dirbpq << std::endl; );
    return compute_line_dir(xpq, dirbpq);
  }

  // check whether the point p is an endpoint of the segment s
  inline static
  bool is_endpoint_of(const Site_2& p, const Site_2& s)
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    Are_same_points_2 same_points;
    return ( same_points(p, s.source_site()) ||
	     same_points(p, s.target_site())   );
  }

  // Orient the segment s and return result as a line.
  // Site p is a point which is an endpoint of s and
  // p_before_s is true if p is just before s in the
  // Voronoi vertex.
  static Line_2
  orient_line_endp(const Site_2& p, const Site_2& s, const bool p_before_s)
  {
    return compute_line_from_to(
        p_before_s ? p.point() : other_site(p, s).point(),
        p_before_s ? other_site(p, s).point() : p.point()  );
  }

  // Orient the segment s and return result as a line.
  // Site p is a point which is not an endpoint of s.
  // Site p must not be on the line defined by s.
  static Line_2
  orient_line_nonendp(const Site_2& p, const Site_2& s)
  {
    Line_2 lseg = compute_supporting_line(s.supporting_site());
    Oriented_side os = oriented_side_of_line(lseg, p.point());
    if (os != ON_POSITIVE_SIDE) {
      CGAL_assertion( os == ON_NEGATIVE_SIDE );
      lseg = opposite_line(lseg);
    }
    return lseg;
  }

  // given that point site p is an endpoint of segment seg,
  // return the other endpoint of seg (as a site)
  inline static
  const Site_2 other_site(const Site_2& p, const Site_2& seg)
  {
    CGAL_precondition( p.is_point() && seg.is_segment() );
    Are_same_points_2 same_points;
    if ( same_points(p, seg.source_site()) ){
      return seg.target_site();
    } else {
      CGAL_assertion(same_points(p, seg.target_site()));
      return seg.source_site();
    }
  }

  // Given corner with bearing cb (0: bottom right, 2: top right,
  // 4: top left, 6: bottom left) and another point p, return the
  // (smallest) square having that corner and passing through p.
  static
  Point_2 center_from_corner_and_pt(
      const Point_2 & corner, const Bearing cb, const Point_2 & p)
  {
    CGAL_precondition(cb % 2 == 0);
    const FT absdifx = CGAL::abs(corner.x() - p.x());
    const FT absdify = CGAL::abs(corner.y() - p.y());
    const Comparison_result cmp = CGAL::compare(absdifx, absdify);
    if (cmp == SMALLER) {
      const FT ox = corner.x() + FT((cb < 3) ? -1: +1)*absdify/FT(2);
      const FT oy = (corner.y() + p.y())/FT(2);
      return Point_2(ox, oy);
    } else {
      const FT ox = (corner.x() + p.x())/FT(2);
      const FT oy = corner.y() + FT((cb % 6 == 0) ? +1: -1)*absdifx/FT(2);
      return Point_2(ox, oy);
    }
  }

  inline
  static
  Point_2 center_from_opposite_corners(
      const Point_2 & c, const Point_2 & d)
  {
    return Point_2(c.x() + d.x(), c.y() + d.y(), RT(2));
  }

  inline
  static
  Point_2 center_from_same_side_corners(
      const Point_2 & c, const Point_2 & d, const Bearing bside)
  {
    CGAL_precondition(bside % 2 == 1);
    const FT ax = (bside % 4 == 1) ?
      RT(2)*c.x() + c.y() - d.y() : c.x() + d.x();
    const FT ay = (bside % 4 == 1) ?
      c.y() + d.y() : RT(2)*c.y() + d.x() - c.y();
    return Point_2(ax, ay, RT(2));
  }

  inline
  static
  bool points_inside_touching_sides_v(
      const Line_2 & ls, const Site_2 & pt_site,
      const Site_2 & other_s, const Site_2 & t, const Point_2 & v)
  {
    CGAL_precondition(pt_site.is_point());
    CGAL_precondition(t.is_point());
    CGAL_USE(other_s);
    CGAL_SDG_DEBUG(std::cout << "debug points_inside_touching_sides_v "
        << "ls: " << ls.a() << ' ' << ls.b() << ' ' <<  ls.c()
        << " pt_site=" << pt_site << " other_s=" << other_s
        << " t=" << t << " v=" << v << std::endl;);
    const Point_2 corner =
      compute_linf_projection_nonhom(ls, v);
    const Line_2 ltest = has_positive_slope(ls) ?
      compute_pos_45_line_at(v): compute_neg_45_line_at(v);
    CGAL_assertion(
        oriented_side_of_line(ltest, v) == ON_ORIENTED_BOUNDARY);
    const Oriented_side ost = oriented_side_of_line(ltest, t.point());
    const Oriented_side osx = oriented_side_of_line(ltest, corner);
    CGAL_SDG_DEBUG(std::cout << "debug points_inside_touching_sides_v"
        << " ltest: " << ltest.a() << ' ' << ltest.b() << ' ' <<  ltest.c()
        << " v=" << v << " ost=" << ost
        << " corner=" << corner << " osx=" << osx << std::endl;);
    if (ost == osx) {
      const Point_2 & p = pt_site.point();
      const Oriented_side osp = oriented_side_of_line(ltest, p);
      if (ost == osp) {
        // +-pi/2 slope line through corner and v
        const Line_2 lcv = has_positive_slope(ls) ?
          compute_neg_45_line_at(v): compute_pos_45_line_at(v);
        const Oriented_side oslt = oriented_side_of_line(lcv, t.point());
        const Oriented_side oslp = oriented_side_of_line(lcv, p);
        if (oslt != oslp) {
          return true;
        }
      }
    }
    return false;
  }

  inline
  static
  bool points_inside_touching_sides_v(
      const Site_2 & s, const Site_2 & pt_site,
      const Site_2 & other_s, const Site_2 & t, const Point_2 & v)
  {
    CGAL_precondition(t.is_point());
    CGAL_precondition(pt_site.is_point());
    CGAL_precondition(s.is_segment());
    CGAL_SDG_DEBUG(std::cout << "debug points_inside_touching_sides_v "
        << "s=" << s
        << " pt_site=" << pt_site << " other_s=" << other_s
        << " t=" << t << std::endl;);
    CGAL_assertion(! is_site_h_or_v(s));
    if (other_s.is_segment()) {
      // shortcut: when the point pt_site is on a corner of
      // the Linf square, because it is the endpoint of the
      // other site which is a segment; return false immediately
      if ((! is_site_h_or_v(other_s)) &&
          is_endpoint_of(pt_site, other_s)) {
        return false;
      }
    }
    const Line_2 ls = compute_supporting_line(s.supporting_site());
    return points_inside_touching_sides_v(ls, pt_site, other_s, t, v);
  }

  // Check if point p's Voronoi area has zero area, because point p
  // is sandwiched between two segments with agreeing slope and
  // direction. As a result: vpqr and vqps coincide in edge conflicts.
  inline static
  bool
  zero_voronoi_area(const Site_2& p, const Site_2& r, const Site_2& s)
  {
    Are_same_points_2 same_points;
    if (p.is_segment()) { return false; }
    if (r.is_point() || s.is_point()) { return false; }
    const bool is_p_rsrc = same_points(p, r.source_site());
    const bool is_p_rtrg =
      (! is_p_rsrc) && same_points(p, r.target_site());
    const bool is_p_endp_of_r = is_p_rsrc || is_p_rtrg;
    if (is_p_endp_of_r) {
      const bool is_p_ssrc = same_points(p, s.source_site());
      const bool is_p_strg =
        (! is_p_ssrc) && same_points(p, s.target_site());
      const bool is_p_endp_of_s = is_p_ssrc || is_p_strg;
      if (is_p_endp_of_s) {
        if (is_site_horizontal(r) && is_site_horizontal(s)) { return true; }
        if (is_site_vertical(r) && is_site_vertical(s)) { return true; }
        if ((! is_site_h_or_v(r)) && (! is_site_h_or_v(s))) {
          const bool pos_r = has_positive_slope(r);
          const bool pos_s = has_positive_slope(s);
          if (pos_r == pos_s) {
            const Line_2 l = pos_r ? compute_neg_45_line_at(p.point()) :
              compute_pos_45_line_at(p.point()) ;
            const Oriented_side osr =
              oriented_side_of_line(l, is_p_rsrc ? r.target() : r.source());
            const Oriented_side oss =
              oriented_side_of_line(l, is_p_ssrc ? s.target() : s.source());
            if (osr != oss) { return true; }
          }
        }
      }
    }
    return false;
  }

  // Check if point p is on the line of the axis-parallel segment s.
  // It returns true only for an axis-parallel segment s argument.
  static inline bool is_on_hv_seg_line(const Site_2 & p, const Site_2 & s)
  {
    CGAL_precondition(p.is_point());
    CGAL_precondition(s.is_segment());
    Compare_x_2_Sites_Type scmpx;
    Compare_y_2_Sites_Type scmpy;
    const bool is_hor = is_site_horizontal(s);
    const bool is_ver = (! is_hor) && is_site_vertical(s);
    if (is_hor || is_ver) {
      return ( (is_hor) ?
               scmpy(p, s.source_site()) : scmpx(p, s.source_site()) )
          == EQUAL;
    } else {
      return false;
    }
  }


}; // end of struct Basic_predicates_C2


} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_BASIC_PREDICATES_C2_H

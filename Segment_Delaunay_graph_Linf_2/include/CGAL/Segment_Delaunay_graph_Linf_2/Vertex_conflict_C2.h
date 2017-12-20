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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VERTEX_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VERTEX_CONFLICT_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Orientation_Linf_2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//---------------------------------------------------------------------

template<class K, class Method_tag>
class Vertex_conflict_C2
  : public Basic_predicates_C2<K>
{
private:

  typedef Basic_predicates_C2<K> Base;
  using Base::compute_supporting_line;
  using Base::oriented_side_of_line;
  //using Base::opposite_line;
  using Base::compute_linf_perpendicular;
  using Base::compute_line_from_to;
  using Base::compute_horizontal_line_from_to;
  using Base::compute_vertical_line_from_to;
  using Base::compute_horizontal_side_line;
  using Base::compute_vertical_side_line;
  using Base::compute_perpendicular;
  using Base::compute_cw_perpendicular;
  using Base::compute_horizontal_projection;
  using Base::compute_vertical_projection;
  using Base::intersects_segment_interior_inf_box;
  using Base::intersects_segment_positive_of_wedge;
  using Base::intersects_segment_negative_of_wedge;
  using Base::has_positive_slope;
  using Base::is_site_h_or_v;

  typedef typename Base::Line_2              Line_2;

  typedef typename K::Point_2                Point_2;
  typedef typename K::Segment_2              Segment_2;
  typedef typename K::Site_2                 Site_2;
  typedef typename K::FT                     FT;
  typedef typename K::RT                     RT;
  //typedef typename K::Orientation            Orientation;
  typedef typename K::Sign                   Sign;

  typedef Voronoi_vertex_C2<K,Method_tag>    Voronoi_vertex_2;

  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>   Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K> Are_same_segments_2;

  typedef typename K::Intersections_tag      ITag;

  // the orientation Linf predicate for three points
  typedef Orientation_Linf_2<K>              Orientation_Linf_points_2;

  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;

  Compare_x_2 cmpx;
  Compare_y_2 cmpy;
  Orientation_Linf_points_2 orientation_linf;


private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

  bool is_on_common_support(const Site_2& s1, const Site_2& s2,
			    const Point_2& p) const
  {
    CGAL_precondition( !s1.is_input() && !s2.is_input() );

    if (  same_segments(s1.supporting_site(0),
			s2.supporting_site(0)) ||
	  same_segments(s1.supporting_site(0),
			s2.supporting_site(1))  ) {
      Site_2 support = s1.supporting_site(0);
      Site_2 tp = Site_2::construct_site_2(p);

      return (  same_points(support.source_site(), tp) ||
		same_points(support.target_site(), tp)  );
    } else if (  same_segments(s1.supporting_site(1),
			       s2.supporting_site(1)) ||
		 same_segments(s1.supporting_site(1),
			       s2.supporting_site(0))  ) {
      Site_2 support = s1.supporting_site(1);
      Site_2 tp = Site_2::construct_site_2(p);

      return (  same_points(support.source_site(), tp) ||
		same_points(support.target_site(), tp)  );
    }
    return false;
  }

  bool have_common_support(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( !p.is_input() && !q.is_input() );

    return
      same_segments(p.supporting_site(0), q.supporting_site(0)) ||
      same_segments(p.supporting_site(0), q.supporting_site(1)) ||
      same_segments(p.supporting_site(1), q.supporting_site(1)) ||
      same_segments(p.supporting_site(1), q.supporting_site(0));
  }

  bool have_common_support(const Site_2& s, const Point_2& p1,
			   const Point_2& p2) const
  {
    CGAL_precondition( !s.is_input() );

    Site_2 t = Site_2::construct_site_2(p1, p2);

    return ( same_segments(s.supporting_site(0), t) ||
	     same_segments(s.supporting_site(1), t) );
  }

private:
  Sign incircle_ppp(const Site_2& p, const Site_2& q,
		    const Site_2& t, const Tag_false&) const
  {
    Point_2 pp = p.point(), qp = q.point(), tp = t.point();

    // MK::ERROR: here I should call a kernel object, not a
    // function...; actually here (and everywhere in this class)
    // use the orientation predicate for sites; it does some
    // geometric filtering...
    // philaris: I think we call an object now
    OrientationLinf o = orientation_linf(pp, qp, tp);

    if ( o != DEGENERATE ) {
      return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
    }

    // here orientation Linf is degenerate

    // return NEGATIVE if t is strictly inside the bounding box
    // of p and q, otherwise return POSITIVE;
    // ZERO is not possible (because t should not be equal
    // to either p or q)

    Comparison_result cmpxpt = cmpx(pp, tp);
    Comparison_result cmpxtq = cmpx(tp, qp);
    Comparison_result cmpypt = cmpy(pp, tp);
    Comparison_result cmpytq = cmpy(tp, qp);

    //Sign s = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
    Sign s = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

    CGAL_assertion( s != ZERO );

    return s;
  }

  Sign incircle_ppp(const Site_2& p, const Site_2& q,
		    const Site_2& t, const Tag_true&) const
  {
    Orientation o = DEGENERATE; // the initialization was done in
                               // order a compiler warning

    // do some geometric filtering...
    bool p_exact = p.is_input();
    bool q_exact = q.is_input();
    bool t_exact = t.is_input();
    bool filtered = false;
    // the following if-statement does the gometric filtering...
    // maybe it is not so important since this will only be
    // activated if a lot of intersection points appear on the
    // convex hull
    if ( !p_exact || !q_exact || !t_exact ) {
      if ( !p_exact && !q_exact && !t_exact ) {
	if ( have_common_support(p, q) &&
	     have_common_support(q, t) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( !p_exact && !q_exact && t_exact ) {
	if ( is_on_common_support(p, q, t.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( !p_exact && q_exact && !t_exact ) {
	if ( is_on_common_support(p, t, q.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( p_exact && !q_exact && !t_exact ) {
	if ( is_on_common_support(t, q, p.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( !p_exact && q_exact && t_exact ) {
	if ( have_common_support(p, q.point(), t.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( p_exact && !q_exact && t_exact ) {
	if ( have_common_support(q, p.point(), t.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( p_exact && q_exact && !t_exact ) {
	if ( have_common_support(t, p.point(), q.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      }
    }

    Point_2 pp = p.point(), qp = q.point(), tp = t.point();

    if ( !filtered ) {
      // MK::ERROR: here I should call a kernel object, not a
      // function...; actually here (and everywhere in this class)
      // use the orientation predicate for sites; it does some
      // geometric filtering...
      // philaris: I think we call an object now
      o = orientation_linf(pp, qp, tp);
    }

    if ( o != DEGENERATE ) {
      return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
    }

    // here orientation Linf is degenerate

    // return NEGATIVE if t is strictly inside the bounding box
    // of p and q, otherwise return POSITIVE;
    // ZERO is not possible (because t should not be equal
    // to either p or q)

    Comparison_result cmpxpt = cmpx(pp, tp);
    Comparison_result cmpxtq = cmpx(tp, qp);
    Comparison_result cmpypt = cmpy(pp, tp);
    Comparison_result cmpytq = cmpy(tp, qp);

    //Sign s = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
    Sign s = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

    CGAL_assertion( s != ZERO );

    return s;
  }


  // check if (p, inf, q) conflicts with t
  Sign incircle_p(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_point() );

    if ( p.is_point() && q.is_point() ) {

#if 1
      return incircle_ppp(p, q, t, ITag());

#else
      Orientation o = DEGENERATE; // the initialization was done in
                                 // order a compiler warning

      // do some geometric filtering...
      bool p_exact = p.is_input();
      bool q_exact = q.is_input();
      bool t_exact = t.is_input();
      bool filtered = false;
      // the following if-statement does the gometric filtering...
      // maybe it is not so important since this will only be
      // activated if a lot of intersection points appear on the
      // convex hull
      if ( !p_exact || !q_exact || !t_exact ) {
	if ( !p_exact && !q_exact && !t_exact ) {
	  if ( have_common_support(p, q) &&
	       have_common_support(q, t) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( !p_exact && !q_exact && t_exact ) {
	  if ( is_on_common_support(p, q, t.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( !p_exact && q_exact && !t_exact ) {
	  if ( is_on_common_support(p, t, q.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( p_exact && !q_exact && !t_exact ) {
	  if ( is_on_common_support(t, q, p.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( !p_exact && q_exact && t_exact ) {
	  if ( have_common_support(p, q.point(), t.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( p_exact && !q_exact && t_exact ) {
	  if ( have_common_support(q, p.point(), t.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( p_exact && q_exact && !t_exact ) {
	  if ( have_common_support(t, p.point(), q.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	}
      }

      Point_2 pp = p.point(), qp = q.point(), tp = t.point();

      if ( !filtered ) {
	// MK::ERROR: here I should call a kernel object, not a
	// function...; actually here (and everywhere in this class)
	// use the orientation predicate for sites; it does some
	// geometric filtering...
	// philaris: I think we call an object now
	o = orientation_linf(pp, qp, tp);
      }

      if ( o != DEGENERATE ) {
	return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
      }

      // here orientation Linf is degenerate

      // return NEGATIVE if t is strictly inside the bounding box
      // of p and q, otherwise return POSITIVE;
      // ZERO is not possible (because t should not be equal
      // to either p or q)

      Comparison_result cmpxpt = cmpx(pp, tp);
      Comparison_result cmpxtq = cmpx(tp, qp);
      Comparison_result cmpypt = cmpy(pp, tp);
      Comparison_result cmpytq = cmpy(tp, qp);

      //Sign s = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
      Sign s = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

      CGAL_assertion( s != ZERO );

      return s;
#endif
    }

    // here the pair of p, q contains a point and a segment

    CGAL_assertion( p.is_point() || q.is_point() );

    CGAL_SDG_DEBUG( std::cout << "debug incircle_p of p=" << p
                   << " q=" << q << " t=" << t << std::endl; );

    // philaris: here, for Linf there are differences from L2

    if ((q.is_segment() && is_site_h_or_v(q)) ||
        (p.is_segment() && is_site_h_or_v(p))   ) {
      // the case of horizontal or vertical segment is the same
      // as in L2
      Orientation o;
      if ( p.is_point() && q.is_segment() ) {
        CGAL_assertion( same_points(p, q.source_site()) ||
            same_points(p, q.target_site()) );
        Point_2 pq =
          same_points(p, q.source_site()) ? q.target() : q.source();
        o = orientation_linf(p.point(), pq, t.point());
      } else { // p is a segment and q is a point
        CGAL_assertion( same_points(q, p.source_site()) ||
            same_points(q, p.target_site()) );
        Point_2 pp =
          same_points(q, p.source_site()) ? p.target() : p.source();
        o = orientation_linf(pp, q.point(), t.point());
      }
      if ( CGAL::is_certain(o == RIGHT_TURN) )  {
        CGAL_SDG_DEBUG(std::cout << "debug incircle_p about to return "
          << (CGAL::get_certain( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE)
                       << std::endl ;);
        return CGAL::get_certain( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE;
      }
      return CGAL::Uncertain<CGAL::Sign>::indeterminate();
    } else {
      // here the segment is neither horizontal nor vertical
      if ( p.is_point() && q.is_segment() ) {

        CGAL_SDG_DEBUG(std::cout
            << "debug incircle_p: p point, q segment" << std::endl; );

        Point_2 pnt = p.point();

        // compute slope of segment q
        Line_2 l = compute_supporting_line(q.supporting_site());
        bool is_positive_slope =
          CGAL::sign(l.a()) != CGAL::sign(l.b());

        CGAL_SDG_DEBUG(std::cout << "debug incircle_p: q has "
          << (is_positive_slope ? "positive" : "negative")
          << " slope" << std::endl; );


        Oriented_side side_of_pnt;
        if (same_points(q.source_site(), p) ||
            same_points(q.target_site(), p)   ) {
          side_of_pnt = ON_ORIENTED_BOUNDARY;
        } else {
          side_of_pnt = oriented_side_of_line(l, pnt);
        }

        if (side_of_pnt == ON_ORIENTED_BOUNDARY) {
          CGAL_assertion(same_points(q.source_site(), p) ||
                         same_points(q.target_site(), p)   );
        } else {
          // here, side_of_pnt is non-zero
          Oriented_side side_of_t = oriented_side_of_line(l, t.point());
          if (side_of_t == - side_of_pnt) {
            // p and t on different sides of segment q
            return POSITIVE;
          }
        }

        Line_2 lhor;
        Line_2 lver;

        bool is_same_qsrc_p = same_points(q.source_site(), p);
        bool is_same_qtrg_p = same_points(q.target_site(), p);

        Segment_2 seg = q.segment();

        if (is_same_qsrc_p || is_same_qtrg_p) {
          const Point_2 otherpnt =
            is_same_qsrc_p ? seg.target() : seg.source();
          lhor = compute_horizontal_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                    ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE );
          lver = compute_vertical_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                   ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE);
        } else {
          // here point p is not on segment q
          Point_2 pnt_on_seg;
          if (is_positive_slope) {
            pnt_on_seg = compute_vertical_projection(l, pnt);
            lver = compute_line_from_to(pnt_on_seg, pnt);
            lhor = compute_cw_perpendicular(lver, pnt_on_seg);

            if ( cmpx(pnt, t.point()) == EQUAL ) {
              Comparison_result cytprj = cmpy(t.point(), pnt_on_seg);

              if (cytprj == EQUAL) {
                // here t equals the projection
                return POSITIVE;
              }

              Comparison_result cypntt = cmpy(pnt, t.point());

              if (cypntt == cytprj) {
                return POSITIVE;
              } else {
                return NEGATIVE;
              }
            } // end of case where pnt and t have same x

          } else { // is_negative_slope
            pnt_on_seg = compute_horizontal_projection(l, pnt);
            lhor = compute_line_from_to(pnt_on_seg, pnt);
            lver = compute_cw_perpendicular(lhor, pnt_on_seg);

            if ( cmpy(pnt, t.point()) == EQUAL ) {
              Comparison_result cxtprj = cmpx(t.point(), pnt_on_seg);

              if (cxtprj == EQUAL) {
                // here t equals the projection
                return POSITIVE;
              }

              Comparison_result cxpntt = cmpx(pnt, t.point());

              if (cxpntt == cxtprj) {
                return POSITIVE;
              } else {
                return NEGATIVE;
              }
            } // end of case where pnt and t have same y

          }
        } // end of else of if (is_same_qsrc_p or is_same_qtrg_p)

        // here, use lhor and lver to decide about t
        // philaris: negative means conflict
        //           positive means no conflict

        CGAL_SDG_DEBUG(std::cout
            << "CGAL_SDG_DEBUG: debug incircle_p lhor=("
            << lhor.a() << " " << lhor.b() << " " << lhor.c()
            << ") lver=("
            << lver.a() << " " << lver.b() << " " << lver.c()
                       << ")" << std::endl ;);

        Oriented_side osh =
          oriented_side_of_line(lhor, t.point());

        Oriented_side osv =
          oriented_side_of_line(lver, t.point());

        CGAL_SDG_DEBUG(std::cout << "debug incircle_p osh=" << osh
                       << " osv=" << osv << std::endl; );

        if ((osh !=
             ON_NEGATIVE_SIDE) &&
            (osv !=
             ON_NEGATIVE_SIDE))
        {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_p about to return NEG" << std::endl; );
          return NEGATIVE;
        } else {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_p about to return POS" << std::endl; );
          return POSITIVE;
        }

      } else { // p is a segment and q is a point

        CGAL_SDG_DEBUG(std::cout
            << "debug incircle_p: p segment, q point" << std::endl; );

        // compute slope of segment p
        Line_2 l = compute_supporting_line(p.supporting_site());
        bool is_positive_slope =
          CGAL::sign(l.a()) != CGAL::sign(l.b());

        Point_2 pnt = q.point();

        Oriented_side side_of_pnt;
        if (same_points(p.source_site(), q) ||
            same_points(p.target_site(), q)   ) {
          side_of_pnt = ON_ORIENTED_BOUNDARY;
        } else {
          side_of_pnt = oriented_side_of_line(l, pnt);
        }

        CGAL_SDG_DEBUG(std::cout << "debug incircle_p"
           << " side_of_pnt=" << side_of_pnt << std::endl; );

        if (side_of_pnt == ON_ORIENTED_BOUNDARY) {
          CGAL_assertion(same_points(p.source_site(), q) ||
                         same_points(p.target_site(), q)   );
        } else {
          // here, side_of_pnt is non-zero
          Oriented_side side_of_t = oriented_side_of_line(l, t.point());
          if (side_of_t == - side_of_pnt) {
            // q and t on different sides of segment p
            return POSITIVE;
          }
        }

        Line_2 lhor;
        Line_2 lver;

        bool is_same_psrc_q = same_points(p.source_site(), q);
        bool is_same_ptrg_q = same_points(p.target_site(), q);

        Segment_2 seg = p.segment();

        if (is_same_psrc_q || is_same_ptrg_q) {
          const Point_2 otherpnt = is_same_psrc_q ? seg.target() : seg.source();
          lhor = compute_horizontal_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                    ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE );
          lver = compute_vertical_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                    ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE );
        } else {
          // here point q is not on segment p

          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_p q not on p" << std::endl; );

          Point_2 pnt_on_seg;
          if (is_positive_slope) {
            pnt_on_seg = compute_horizontal_projection(l, pnt);
            lhor = compute_horizontal_line_from_to(pnt, pnt_on_seg);
            lver = compute_perpendicular(lhor, pnt_on_seg);

            if ( cmpy(pnt, t.point()) == EQUAL ) {
              Comparison_result cxtprj = cmpx(t.point(), pnt_on_seg);

              if (cxtprj == EQUAL) {
                // here t equals the projection
                return POSITIVE;
              }

              Comparison_result cxpntt = cmpx(pnt, t.point());

              if (cxpntt == cxtprj) {
                return POSITIVE;
              } else {
                return NEGATIVE;
              }
            } // end of case where pnt and t have same y

          } else { // is_negative_slope
            pnt_on_seg = compute_vertical_projection(l, pnt);
            lver = compute_vertical_line_from_to(pnt, pnt_on_seg);
            lhor = compute_perpendicular(lver, pnt_on_seg);

            if ( cmpx(pnt, t.point()) == EQUAL ) {
              Comparison_result cytprj = cmpy(t.point(), pnt_on_seg);

              if (cytprj == EQUAL) {
                // here t equals the projection
                return POSITIVE;
              }

              Comparison_result cypntt = cmpy(pnt, t.point());

              if (cypntt == cytprj) {
                return POSITIVE;
              } else {
                return NEGATIVE;
              }
            } // end of case where pnt and t have same x

          }
        } // end of else of if (is_same_psrc_q or is_same_ptrg_q)

        // here, use lhor and lver to decide about t
        // philaris: negative means conflict
        //           positive means no conflict

        CGAL_SDG_DEBUG(std::cout
            << "CGAL_SDG_DEBUG: debug vc incircle_p lhor=("
            << lhor.a() << " " << lhor.b() << " " << lhor.c()
            << ") lver=("
            << lver.a() << " " << lver.b() << " " << lver.c()
                       << ")" << std::endl ;);

        Oriented_side osh =
          oriented_side_of_line(lhor, t.point());

        Oriented_side osv =
          oriented_side_of_line(lver, t.point());

        CGAL_SDG_DEBUG(std::cout << "debug vc incircle_p osh=" << osh
            << " osv=" << osv << std::endl; );

        if ((osh !=
             ON_NEGATIVE_SIDE) &&
            (osv !=
             ON_NEGATIVE_SIDE))
        {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_p about to return NEG"
              << std::endl; );
          return NEGATIVE;
        } else {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_p about to return POS"
              << std::endl; );
          return POSITIVE;
        }

      } // end of case where p is a segment and q is a point
    } // end of case where segment is neither horizontal nor vertical
  }

  //-----------------------------------------------------------------------

  // philaris: here Linf is very different from L2
  Sign incircle_pps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_pps p=" << p
      << " q=" << q << " t=" << t << std::endl; );

    const bool is_p_tsrc = same_points(p, t.source_site());
    const bool is_p_ttrg = same_points(p, t.target_site());

    const bool is_q_tsrc = same_points(q, t.source_site());
    const bool is_q_ttrg = same_points(q, t.target_site());

    const bool is_p_on_t = is_p_tsrc || is_p_ttrg;
    const bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    if ( is_p_on_t && is_q_on_t ) {
	// if t is the segment joining p and q then t must be a vertex
	// on the convex hull
	return NEGATIVE;
    } else {
      // here, neither p nor q is on t
      if (is_site_h_or_v(t)) {
        // here segment is horizontal or vertical;
        // there are similarities with L2 case here

        if ( is_p_on_t ) {
          // p is an endpoint of t
          // in this case the p,q,oo vertex is destroyed only if the
          // other endpoint of t is beyond
          Point_2 pt = is_p_tsrc ? t.target() : t.source();
          Orientation o = CGAL::orientation(p.point(), q.point(), pt);
          CGAL_SDG_DEBUG(std::cout << "debug incircle_pps p is endp of t: "
                         << "p = " << p << ", q = " << q << ", t = " << t
                         << " o = " << o << std::endl; );
          return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
        } else if ( is_q_on_t ) {
          Point_2 pt = is_q_tsrc ? t.target() : t.source();
          Orientation o = CGAL::orientation(p.point(), q.point(), pt);
          CGAL_SDG_DEBUG(std::cout << "debug incircle_pps q is endp of t: "
                         << "p = " << p << ", q = " << q << ", t = " << t
                         << " o = " << o << std::endl; );
          return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
        } else {
          // We immediately return POSITIVE here.
          // This only works because the segment t is hor. or ver. !
          // If t is neither hor. or ver., we have a more complicated
          // decision procedure (see below).
          // Since we insert endpoints of segments first, p and q cannot
          // be consecutive points on the convex hull, because of the
          // endpoints of t.
          return POSITIVE;
        }
      } else {
        // here segment is neither horizontal nor vertical

        // consider corner of infinite square going through
        // p and q (if you traverse this square ccw, you meet
        // q first, then p)

        Point_2 pp = p.point();
        Point_2 qq = q.point();

        Comparison_result cmpxpq = cmpx(pp, qq);
        Comparison_result cmpypq = cmpy(pp, qq);

        if ((cmpxpq == EQUAL) || (cmpypq == EQUAL)) {
          if (! (is_p_tsrc || is_q_tsrc)) {
            if (CGAL::right_turn(pp, qq, t.segment().source())) {
              return NEGATIVE;
            }
          }
          if (! (is_p_ttrg || is_q_ttrg)) {
            if (CGAL::right_turn(pp, qq, t.segment().target())) {
              return NEGATIVE;
            }
          }
          return POSITIVE;
        } else {
          // here, p and q do not have any same coordinate
          if (intersects_segment_interior_inf_box(t, q, p, cmpxpq, cmpypq)) {
            return NEGATIVE;
          } else {
            return POSITIVE;
          }
        } // end of case where p and q do not have any same coordinate

        // philaris: execution should never reach here
        // philaris: code added to avoid warnings on some compilers
        CGAL_assertion( false );
        return ZERO;

      } // end of case where segment is neither horizontal nor vertical
    } // end of case where neither p nor q is on t
  } // end of incircle_pps


  Sign incircle_sps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_segment() && q.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_sps p=" << p
      << " q=" << q << " t=" << t << std::endl; );

    bool is_q_tsrc = same_points(q, t.source_site());
    bool is_q_ttrg = same_points(q, t.target_site());

    bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    bool is_q_psrc = same_points(q, p.source_site());
    bool is_q_ptrg = same_points(q, p.target_site());

    bool is_q_on_p = is_q_psrc || is_q_ptrg;


    if ( is_q_on_t && is_q_on_p) {
      CGAL_SDG_DEBUG(std::cout << "debug incircle_sps "
        << "is_q_on_t && is_q_on_p" << std::endl; );

      Point_2 pother = is_q_psrc ? p.target() : p.source();
      Point_2 tother = is_q_tsrc ? t.target() : t.source();

      Orientation o = CGAL::orientation(pother, q.point(), tother);

      CGAL_SDG_DEBUG(std::cout << "debug incircle_sps or( "
        << pother << ", " << q.point() << ", " << tother << " ) = "
        << o << std::endl; );

      if (o != RIGHT_TURN) {
        return POSITIVE;
      } else {
        if (is_site_h_or_v(p) || is_site_h_or_v(t)) {
          return NEGATIVE;
        } else {
          bool has_p_pos_slope = has_positive_slope(p);
          if (has_p_pos_slope) {
            return cmpy(pother, q.point()) == cmpy(tother, q.point()) ?
                   NEGATIVE : POSITIVE;
          } else {
            return cmpx(pother, q.point()) == cmpx(tother, q.point()) ?
                   NEGATIVE : POSITIVE;
          }
        }
      }
    } else {
      // philaris: serious difference for Linf here, related to L2

      if ( is_site_h_or_v(p) ) {
        return POSITIVE;
      } else {
        CGAL_assertion(! is_site_h_or_v(p));

        // compute slope of segment p
        Line_2 l = compute_supporting_line(p.supporting_site());
        bool is_positive_slope =
          CGAL::sign(l.a()) != CGAL::sign(l.b());

        Point_2 pnt = q.point();

        Oriented_side side_of_pnt =
          oriented_side_of_line(l, pnt);

        if (side_of_pnt == ON_ORIENTED_BOUNDARY) {
          CGAL_assertion(same_points(p.source_site(), q) ||
                         same_points(p.target_site(), q)   );
        }

        Line_2 lhor;
        Line_2 lver;

        bool is_same_psrc_q = same_points(p.source_site(), q);
        bool is_same_ptrg_q = same_points(p.target_site(), q);

        Point_2 pnt_on_seg;
        Segment_2 seg = p.segment();

        if (is_same_psrc_q || is_same_ptrg_q) {
          pnt_on_seg = pnt;
          Point_2 otherpnt;
          if (is_same_psrc_q) {
            otherpnt = seg.target();
          } else if (is_same_ptrg_q) {
            otherpnt = seg.source();
          } else {
            // unreachable
            CGAL_assertion( false );
          }

          lhor = compute_horizontal_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                    ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE );
          lver = compute_vertical_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                    ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE );
        } else {
          // here point q is not on segment p

          CGAL_SDG_DEBUG(std::cout
              << "debug: here point q is not on segment p"
              << std::endl; );

          if (is_positive_slope) {
            pnt_on_seg = compute_horizontal_projection(l, pnt);
            lhor = compute_horizontal_line_from_to(pnt, pnt_on_seg);
            lver = compute_perpendicular(lhor, pnt_on_seg);
          } else { // is_negative_slope
            pnt_on_seg = compute_vertical_projection(l, pnt);
            CGAL_SDG_DEBUG(std::cout
                << "debug: pnt_on_seg = " << pnt_on_seg
                << std::endl; );
            lver = compute_vertical_line_from_to(pnt, pnt_on_seg);
            lhor = compute_perpendicular(lver, pnt_on_seg);
          }
        } // end of else of if (is_same_psrc_q or is_same_ptrg_q)

        // here, use lhor and lver to decide about t
        // philaris: negative means conflict
        //           positive means no conflict

        CGAL_SDG_DEBUG(std::cout << "debug: sps "
          << " lhor=" << lhor.a() << " " << lhor.b() << " " << lhor.c()
          << " lver=" << lver.a() << " " << lver.b() << " " << lver.c()
          << std::endl; );

        // philaris: careful here
        //if (intersects_segment_positive_halfplane(t, lhor) and
        //    intersects_segment_positive_halfplane(t, lver)    )
        if (intersects_segment_positive_of_wedge(t, lhor, lver))
        {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_sps about to return NEG"
              << std::endl; );
          return NEGATIVE;
        } else {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_sps does not cross wedge, "
              << "check for common endpoint of p and t"
              << std::endl; );

          bool is_psrc_tsrc =
            same_points(p.source_site(), t.source_site());
          bool is_psrc_ttrg =
            same_points(p.source_site(), t.target_site());
          bool is_ptrg_tsrc =
            same_points(p.target_site(), t.source_site());
          bool is_ptrg_ttrg =
            same_points(p.target_site(), t.target_site());

          Point_2 testp, otherp;

          if (is_psrc_tsrc || is_ptrg_tsrc) {
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_sps tsrc endp of p"
                << std::endl; );

            if (same_points(q, t.target_site())) {
              if ( is_site_h_or_v(t) ) {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_sps horver-non seg comm"
                    << std::endl; );
                return POSITIVE;
              }
            }

            testp = t.source_site().point();
            otherp = t.target_site().point();
          } else if (is_psrc_ttrg || is_ptrg_ttrg) {
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_sps ttrg endp of p"
                << std::endl; );

            if (same_points(q, t.source_site())) {
              if ( is_site_h_or_v(t) ) {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_sps horver-non seg comm"
                    << std::endl; );
                return POSITIVE;
              }
            }

            testp = t.target_site().point();
            otherp = t.source_site().point();
          } else {
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_sps fail endp, return POS"
                << std::endl; );
            return POSITIVE;
          }

          // here p and t have common endpoint testp

          CGAL_SDG_DEBUG(std::cout << "debug testp=" << testp << " pnt="
            << pnt << " pnt_on_seg=" << pnt_on_seg << std::endl; );

          // check if testp equals pnt_on_seg

          if ((cmpx(testp, pnt_on_seg) == EQUAL ) &&
              (cmpy(testp, pnt_on_seg) == EQUAL )   ) {
            // here testp is the same as pnt_on_seg

            Oriented_side oscandidate =
              oriented_side_of_line(is_positive_slope? lhor : lver,
                                    otherp);

            if (oscandidate == ON_POSITIVE_SIDE) {

              Oriented_side osq = oriented_side_of_line(l, q.point());
              Oriented_side ost = oriented_side_of_line(l, otherp);

              CGAL_assertion(osq != ON_ORIENTED_BOUNDARY);

              CGAL_SDG_DEBUG(std::cout << "debug p=" << p << " q=" << q
                  << " t=" << t << " osq="<< osq
                  << " ost=" << ost << std::endl;);

              if (osq == ost) {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_sps sameside return NEG"
                    << std::endl; );
                return NEGATIVE;
              } else {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_sps diffside return POS"
                    << std::endl; );
                return POSITIVE;
              }
            } else {
              CGAL_SDG_DEBUG(std::cout
                  << "debug incircle_sps NOT othert on POS return POS"
                  << std::endl; );
              return POSITIVE;
            }
          } else {
            // here, testp is not the same as pnt_on_seg
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_sps diffpnts return POS"
                << std::endl; );
            return POSITIVE;
          }

        } // end of case where segment does not intersect
          // positive part of wedge

      } // end of case where p is non-hor and non-ver
    } // end of case where q is not on t
  } // end of function incircle_sps


  Sign incircle_pss(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_pss (piqt) = ("
        << p << ") inf (" << q << ") (" << t << ")" << std::endl; );

    bool is_p_tsrc = same_points(p, t.source_site());
    bool is_p_ttrg = same_points(p, t.target_site());

    bool is_p_on_t = is_p_tsrc || is_p_ttrg;

    bool is_p_qsrc = same_points(p, q.source_site());
    bool is_p_qtrg = same_points(p, q.target_site());

    bool is_p_on_q = is_p_qsrc || is_p_qtrg;

    if ( is_p_on_t && is_p_on_q ) {
      CGAL_SDG_DEBUG(std::cout << "debug incircle_pss "
        << "is_p_on_t && is_p_on_q" << std::endl; );

      Point_2 qother = is_p_qsrc ? q.target() : q.source();
      Point_2 tother = is_p_tsrc ? t.target() : t.source();

      Orientation o = CGAL::orientation(qother, p.point(), tother);

      CGAL_SDG_DEBUG(std::cout << "debug incircle_pss or( "
        << qother << ", " << p.point() << ", " << tother << " ) = "
        << o << std::endl; );

      if (o != LEFT_TURN) {
        return POSITIVE;
      } else {
        if ( is_site_h_or_v(q) || is_site_h_or_v(t) ) {
          return NEGATIVE;
        } else {
          bool has_q_pos_slope = has_positive_slope(q);
          if (has_q_pos_slope) {
            return cmpx(qother, p.point()) == cmpx(tother, p.point()) ?
                   NEGATIVE : POSITIVE;
          } else {
            return cmpy(qother, p.point()) == cmpy(tother, p.point()) ?
                   NEGATIVE : POSITIVE;
          }
        }
      }

      return (o == LEFT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // philaris: serious difference for Linf here, related to L2

      if ( is_site_h_or_v(q) ) {
        return POSITIVE;
      } else {
        CGAL_assertion(! is_site_h_or_v(q));

        Point_2 pnt = p.point();

        // compute slope of segment q
        Line_2 l = compute_supporting_line(q.supporting_site());
        bool is_positive_slope =
          CGAL::sign(l.a()) != CGAL::sign(l.b());

        CGAL_SDG_DEBUG(std::cout << "debug incircle_pss: q has " <<
          (is_positive_slope ? "positive" : "negative") <<
          " slope" << std::endl; );

        Oriented_side side_of_pnt =
          oriented_side_of_line(l, pnt);

        if (side_of_pnt == ON_ORIENTED_BOUNDARY) {
          CGAL_assertion(same_points(q.source_site(), p) ||
                         same_points(q.target_site(), p)   );
        }

        Line_2 lhor;
        Line_2 lver;

        bool is_same_qsrc_p = same_points(q.source_site(), p);
        bool is_same_qtrg_p = same_points(q.target_site(), p);

        Point_2 pnt_on_seg;
        Segment_2 seg = q.segment();

        if (is_same_qsrc_p || is_same_qtrg_p) {
          CGAL_SDG_DEBUG(std::cout
              << "debug: pss: p is endpoint of q" << std::endl; );
          pnt_on_seg = pnt;
          Point_2 otherpnt;
          if (is_same_qsrc_p) {
            otherpnt = seg.target();
          } else if (is_same_qtrg_p) {
            otherpnt = seg.source();
          } else {
            // unreachable
            CGAL_assertion( false );
          }

          lhor = compute_horizontal_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                    ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE );
          lver = compute_vertical_side_line(
                   pnt, otherpnt,
                   is_positive_slope ?
                   ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE);
        } else {
          // here point p is not on segment q

          CGAL_SDG_DEBUG(std::cout
              << "debug: pss from pt to pt on seg" << std::endl; );

          if (is_positive_slope) {
            pnt_on_seg = compute_vertical_projection(l, pnt);
            lver = compute_vertical_line_from_to(pnt, pnt_on_seg);
            lhor = compute_cw_perpendicular(lver, pnt_on_seg);
          } else { // is_negative_slope
            pnt_on_seg = compute_horizontal_projection(l, pnt);
            lhor = compute_horizontal_line_from_to(pnt, pnt_on_seg);
            lver = compute_cw_perpendicular(lhor, pnt_on_seg);
          }
        } // end of else of if (is_same_qsrc_p or is_same_qtrg_p)

        // here, use lhor and lver to decide about t
        // philaris: negative means conflict
        //           positive means no conflict

        CGAL_SDG_DEBUG(std::cout << "debug incircle_pss lhor=("
          << lhor.a() << " " << lhor.b() << " " << lhor.c()
          << ") lver=("
          << lver.a() << " " << lver.b() << " " << lver.c()
          << ")" << std::endl; );

        // philaris: careful here
        //if (intersects_segment_negative_halfplane(t, lhor) and
        //    intersects_segment_negative_halfplane(t, lver)    )
        if (intersects_segment_negative_of_wedge(t, lhor, lver))
        {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_pss about to return NEG" << std::endl; );
          return NEGATIVE;
        } else {
          CGAL_SDG_DEBUG(std::cout
              << "debug incircle_pss does not cross wedge, "
              << "check for common endpoint of q and t"
              << std::endl; );

          bool is_qsrc_tsrc =
            same_points(q.source_site(), t.source_site());
          bool is_qsrc_ttrg =
            same_points(q.source_site(), t.target_site());
          bool is_qtrg_tsrc =
            same_points(q.target_site(), t.source_site());
          bool is_qtrg_ttrg =
            same_points(q.target_site(), t.target_site());

          Point_2 testt, othert;

          if (is_qsrc_tsrc || is_qtrg_tsrc) {
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_pss tsrc endp of q" << std::endl; );

            if (same_points(p, t.target_site())) {
              if ( is_site_h_or_v(t) ) {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_sps horver-non seg comm"
                    << std::endl; );
                return POSITIVE;
              }
            }

            testt = t.source_site().point();
            othert = t.target_site().point();
          } else if (is_qsrc_ttrg || is_qtrg_ttrg) {
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_pss ttrg endp of q" << std::endl; );

            if (same_points(p, t.source_site())) {
              if ( is_site_h_or_v(t) ) {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_pss horver-non seg comm"
                    << std::endl; );
                return POSITIVE;
              }
            }

            testt = t.target_site().point();
            othert = t.source_site().point();
          } else {
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_pss fail endp, return POS"
                << std::endl; );
            return POSITIVE;
          }

          // here p and t have common endpoint testt

          CGAL_SDG_DEBUG(std::cout << "debug testt=" << testt
            << " pnt_on_seg=" << pnt_on_seg << std::endl; );

          // check if testp equals pnt_on_seg

          if ((cmpx(testt, pnt_on_seg) == EQUAL ) &&
              (cmpy(testt, pnt_on_seg) == EQUAL )   ) {
            // here testt is the same as pnt_on_seg

            Oriented_side oscandidate =
              oriented_side_of_line(is_positive_slope? lhor : lver,
                                    othert);
            CGAL_SDG_DEBUG(std::cout << "debug oscandidate="
                << oscandidate << std::endl; );

            CGAL_SDG_DEBUG(std::cout << "debug vc other=" << othert
                << " pnt_on_seg=" << pnt_on_seg << std::endl; );

            if (oscandidate == ON_POSITIVE_SIDE) {

              Oriented_side osp = oriented_side_of_line(l, p.point());
              Oriented_side ost = oriented_side_of_line(l, othert);

              CGAL_assertion(osp != ON_ORIENTED_BOUNDARY);

              CGAL_SDG_DEBUG(std::cout << "debug p=" << p << " q=" << q
                  << " t=" << t << " osp="<< osp
                  << " ost=" << ost << std::endl;);

              if (osp == ost) {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_pss sameside return NEG"
                    << std::endl; );
                return NEGATIVE;
              } else {
                CGAL_SDG_DEBUG(std::cout
                    << "debug incircle_pss diffside return POS"
                    << std::endl; );
                return POSITIVE;
              }
            } else {
              CGAL_SDG_DEBUG(std::cout
                  << "debug incircle_pss NOT othert on POS return POS"
                  << std::endl; );
              return POSITIVE;
            }
          } else {
            // here, testt is not the same as pnt_on_seg
            CGAL_SDG_DEBUG(std::cout
                << "debug incircle_pss diffpnts return POS"
                << std::endl; );
            return POSITIVE;
          }

        } // end of case where segment does not intersect
          // negative side of wedge

      } // end of case of non-hor and non-ver segment q
    } // end of case where p is not on t
  } // end of function incircle_pss


  // check if (p, inf, q) conflicts with t
  Sign incircle_s(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_segment() );

    if ( p.is_point() && q.is_point() ) {
      return incircle_pps(p, q, t);
    } else if ( p.is_point() && q.is_segment() ) {
      return incircle_pss(p, q, t);
    } else { // p is a segment and q is a point
      return incircle_sps(p, q, t);
    }
  }


public:
  typedef Site_2      argument_type;
  typedef Sign        result_type;


  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& r, const Site_2& t) const
  {

    CGAL_SDG_DEBUG(std::cout
        << "debug: Vertex_conflict entering (pqrt)= ("
        << p << ") (" << q << ") (" << r << ") (" << t << ")"
        << std::endl; );


#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      int np = 0;
      if ( p.is_point() ) ++np;
      if ( q.is_point() ) ++np;
      if ( r.is_point() ) ++np;
      std::string suffix("-failure-log.cin");
      std::string fname;
      if ( np == 3 ) {
	fname = "ppp";
      } else if ( np == 2 ) {
	fname = "pps";
      } else if ( np == 1 ) {
	fname = "pss";
      } else {
	fname = "sss";
      }
      fname += suffix;
      std::ofstream ofs(fname.c_str(), std::ios_base::app);
      ofs.precision(16);
      ofs << p << std::endl;
      ofs << q << std::endl;
      ofs << r << std::endl;
      ofs << t << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

    Voronoi_vertex_2 v(p, q, r);

    CGAL_SDG_DEBUG(
        Point_2 temppt = v.point();
        std::cout
          << "debug Voronoi_vertex_2 v=" << temppt << std::endl; );

    Sign retval = v.incircle(t);

    CGAL_SDG_DEBUG(std::cout << "debug: Vertex_conflict (pqrt)= ("
        << p << ") (" << q << ") ("  << r << ") (" << t << ") returns "
        << retval << std::endl; );

    return retval;
  }




  // check if (p, q, inf) conflicts with t
  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_SDG_DEBUG(std::cout << "debug: Vertex_conflict (pqit)= ("
        << p << ") (" << q << ") inf (" << t << ")"
        << std::endl; );

#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      std::ofstream ofs("failure-log.cin", std::ios_base::app);
      ofs.precision(16);
      ofs << p << std::endl;
      ofs << q << std::endl;
      ofs << t << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    // philaris: the assertion that the point is endpoint of
    // the segment is not true in Linf (in contrast with L2)
    /*
    if ( p.is_point() && q.is_segment() ) {
      // p must be an endpoint of q
      CGAL_assertion( same_points(p, q.source_site()) ||
		      same_points(p, q.target_site()) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( same_points(p.source_site(), q) ||
		      same_points(p.target_site(), q) );
    }
    */

    Sign retval;

    if ( t.is_point() ) {
      //      return incircle_p(p, q, t);
      retval = incircle_p(q, p, t);
    } else {
      // MK::ERROR: do geometric filtering when orientation is called.
      //    return incircle_s(p, q, t);
      retval = incircle_s(q, p, t);
    }

    CGAL_SDG_DEBUG(std::cout << "debug: Vertex_conflict (pqit)= ("
        << p << ") (" << q << ") inf (" << t << ")"
        << " returns " << retval << std::endl; );
    return retval;
  }


};

//---------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VERTEX_CONFLICT_C2_H

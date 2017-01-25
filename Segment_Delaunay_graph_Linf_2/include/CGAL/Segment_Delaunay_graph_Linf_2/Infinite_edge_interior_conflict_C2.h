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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Orientation_Linf_2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//-----------------------------------------------------------------------------

template<class K, class Method_tag>
class Infinite_edge_interior_conflict_C2
  : public Basic_predicates_C2<K>
{
public:

  typedef Basic_predicates_C2<K>              Base;

  typedef typename K::Site_2           Site_2;
  typedef typename K::Point_2          Point_2;
  typedef typename K::Direction_2      Direction_2;
  typedef typename K::RT               RT;
  typedef typename K::Boolean          Boolean;

  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;

  Compare_x_2 cmpx;
  Compare_y_2 cmpy;

  typedef Boolean                      result_type;
  struct argument_type {};

private:
  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>   Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K> Are_same_segments_2;

  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

  typedef Orientation_Linf_2<K>        Orientation_Linf_2_Type;

  Orientation_Linf_2_Type  or_linf;

  using Base::bounded_side_of_bbox;
  using Base::compute_line_from_to;
  using Base::oriented_side_of_line;
  using Base::compute_supporting_line;
  using Base::compute_vertical_projection;
  using Base::compute_horizontal_projection;
  using Base::compute_linf_projection_nonhom;
  using Base::is_site_h_or_v;
  using Base::zero_voronoi_area;

  typedef typename Base::Line_2        Line_2;

public:
  Boolean   operator()(const Site_2& q, const Site_2& s, const Site_2& r,
		       const Site_2& t, Sign sgn) const
  {

    CGAL_SDG_DEBUG(
        std::cout << "debug infinite-edge-int-cf entering (q,s,r,t,sgn)= "
        << q << ' ' << s << ' ' << r << ' ' << t
        << ' ' << sgn << std::endl;);

    if (sgn == NEGATIVE) {
      if (zero_voronoi_area(q, s, r)) { return true; }
    }

    if ( t.is_segment() ) {
      if (q.is_point() && s.is_point() && r.is_point()) {
        if (sgn == NEGATIVE) {
          CGAL_SDG_DEBUG(std::cout << "debug return tocheck" << std::endl;);
          if (same_points(q, t.source_site()) ||
              same_points(q, t.target_site())   ) {
            // this works because endpoints of a segment are
            // inserted before its interior
            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << false << std::endl;);
            return false;
          } else {

            // here q is point
            if ( ! is_site_h_or_v(t) ) {

              Line_2 lt = compute_supporting_line(t.supporting_site());

              // Linf-project point q to line lt

              Point_2 projq =
                compute_linf_projection_nonhom(lt, q.point());

              Line_2 lq = compute_line_from_to(projq, q.point());

              Oriented_side oss =
                oriented_side_of_line(lq, s.point());

              Oriented_side osr =
                oriented_side_of_line(lq, r.point());

              if ((oss == ON_NEGATIVE_SIDE) &&
                  (osr == ON_POSITIVE_SIDE)    ) {
                CGAL_SDG_DEBUG(
                    std::cout
                    << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
                    << q << ' ' << s << ' ' << r << ' ' << t
                    << ' ' << sgn << " returns "
                    << false << std::endl;);
                return false;
              }
            } // case where t is neither horizontal nor vertical

            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << true << std::endl;);
            return true;
          }
        }
      } // end of case where q, s, r are all points

      CGAL_SDG_DEBUG(std::cout
          << "debug tocheck q,s,r not all points" << std::endl;);

      if (sgn == NEGATIVE) {
        CGAL_SDG_DEBUG(std::cout
            << "debug not all pts return true tocheck" << std::endl;);

        if (q.is_point()) {
          bool is_q_tsrc = same_points(q, t.source_site());
          bool is_q_ttrg = same_points(q, t.target_site());

          CGAL_assertion(! (is_q_tsrc || is_q_ttrg));

          if (is_q_tsrc || is_q_ttrg) {
            // philaris: this code should never be executed

            CGAL_assertion_code(
            bool is_q_ssrc = same_points(q, s.source_site());
            bool is_q_strg = same_points(q, s.target_site());
            bool is_q_rsrc = same_points(q, r.source_site());
            bool is_q_rtrg = same_points(q, r.target_site());
            );

            CGAL_assertion(
                (is_q_ssrc || is_q_strg) &&
                (is_q_rsrc || is_q_rtrg)    );

            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf tocheck com (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << true << std::endl;);
            return true;
          }

          CGAL_assertion( ! (is_q_tsrc || is_q_ttrg) );

          // here q is point
          if ( ! is_site_h_or_v(t) ) {

            Line_2 lt = compute_supporting_line(t.supporting_site());

            // Linf-project point q to line lt

            Point_2 projq =
              compute_linf_projection_nonhom(lt, q.point());

            Line_2 lq = compute_line_from_to(projq, q.point());

            Point_2 srep;
            if (s.is_point()) {
              srep = s.point();
            } else {
              // s is segment
              CGAL_assertion( ! is_site_h_or_v(s) ) ;

              Direction_2 d (s.supporting_site().segment());
              Line_2 ls = compute_supporting_line(s.supporting_site());
              if (CGAL::sign(d.dx()) == CGAL::sign(d.dy())) {
                srep = compute_horizontal_projection(
                    ls, q.point());
              } else {
                srep = compute_vertical_projection(
                    ls, q.point());
              }
            }
            Oriented_side oss =
              oriented_side_of_line(lq, srep);

            Point_2 rrep;
            if (r.is_point()) {
              rrep = r.point();
            } else {
              // r is segment
              CGAL_assertion( ! is_site_h_or_v(r) ) ;

              Direction_2 d (r.supporting_site().segment());
              Line_2 lr = compute_supporting_line(r.supporting_site());
              if (CGAL::sign(d.dx()) == CGAL::sign(d.dy())) {
                rrep = compute_vertical_projection(
                    lr, q.point());
              } else {
                rrep = compute_horizontal_projection(
                    lr, q.point());
              }

            }
            Oriented_side osr =
              oriented_side_of_line(lq, rrep);

            if ((oss == ON_NEGATIVE_SIDE) &&
                (osr == ON_POSITIVE_SIDE)    ) {
              CGAL_SDG_DEBUG(
                  std::cout
                  << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
                  << q << ' ' << s << ' ' << r << ' ' << t
                  << ' ' << sgn << " returns "
                  << false << std::endl;);
              return false;
            }
          } // case where t is neither horizontal nor vertical



          //CGAL_assertion(false);

          // philaris: tocheck more

          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << true << std::endl;);

          return true;

        }

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << true << std::endl;);

        return true;
      }

      // philaris: tocheck
      CGAL_SDG_DEBUG(
          std::cout
          << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
          << q << ' ' << s << ' ' << r << ' ' << t
          << ' ' << sgn << " returns "
          << false << std::endl;);
      return false;
    } // end of case where t is segment

    // here and below, t is always a point

    if ( q.is_segment() ) {

      // philaris: difference from L2 here;
      // in L2, r and s are endpoints of q
      // in Linf they still have to be points, but
      // they do not have to be endpoints of q
      // (this has to be checked)
      CGAL_assertion(s.is_point() && r.is_point());

      if ( is_site_h_or_v(q) )
      {
        // in this case r and s must be endpoints of q
        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << ( sgn == NEGATIVE ) << std::endl;);
        return ( sgn == NEGATIVE );
      } 
      else
      {
        if (sgn == NEGATIVE) {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << true << std::endl;);
          return true;
        }

        // here sgn is non-negative

        bool is_s_endp_of_q =
          same_points(s, q.source_site()) ||
          same_points(s, q.target_site());
        bool is_r_endp_of_q =
          same_points(r, q.source_site()) ||
          same_points(r, q.target_site());

        Line_2 l;
        bool is_conflicting_side_of_q = false;

        if (is_s_endp_of_q && is_r_endp_of_q) {
          // check if t is on positive side of rs
          l = compute_line_from_to(r.point(), s.point());
          if (oriented_side_of_line(l, t.point()) ==
              ON_POSITIVE_SIDE) {
            is_conflicting_side_of_q = true;
          }
        } else {
          l = compute_supporting_line(q.supporting_site());
          Oriented_side sidelt =
            oriented_side_of_line(l, t.point());
          if (is_s_endp_of_q) {
            // here r is not endpoint of q
            Oriented_side sidelr =
              oriented_side_of_line(l, r.point());
            CGAL_assertion(sidelr != ON_ORIENTED_BOUNDARY);
            if (sidelt == sidelr) {
              is_conflicting_side_of_q = true;
            }
          } else {
            // here s is not endpoint of q
            Oriented_side sidels =
              oriented_side_of_line(l, s.point());
            CGAL_assertion(sidels != ON_ORIENTED_BOUNDARY);
            if (sidelt == sidels) {
              is_conflicting_side_of_q = true;
            }
          }
        } // end of case: some of s, r is not endpoint of q

        if (! is_conflicting_side_of_q) {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << false << std::endl;);
          return false;
        }

        // here is_conflicting_side_of_q is true;

        CGAL_SDG_DEBUG(std::cout <<
            "debug infcf is_cf_side_of_q" << std::endl;);

        // compute infinite square such that:
        // if you traverse it ccw, then it meets r and then s

        Comparison_result A = cmpx(s.point(), r.point());
        Comparison_result B = cmpy(s.point(), r.point());

        CGAL_assertion((A != EQUAL) && (B != EQUAL));

        // corner point of infinite square
        Point_2 corner;
        if (A == B) {
          corner = Point_2(s.point().x(), r.point().y());
        } else {
          corner = Point_2(r.point().x(), s.point().y());
        }

        Line_2 lr = compute_line_from_to(r.point(), corner);
        Line_2 ls = compute_line_from_to(corner, s.point());

        Oriented_side sidelrt =
          oriented_side_of_line(lr, t.point());

        Oriented_side sidelst =
          oriented_side_of_line(ls, t.point());

        if ((sidelrt == ON_NEGATIVE_SIDE) ||
            (sidelst == ON_NEGATIVE_SIDE)   ) {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << false << std::endl;);
          return false;
        } else {
          // both sidelrt and sidelst are non-negative
          if (sidelrt == ON_ORIENTED_BOUNDARY) {
            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << false << std::endl;);
            return false;
          } else if (sidelst == ON_ORIENTED_BOUNDARY) {
            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << false << std::endl;);
            return false;
          } else {
            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << true << std::endl;);
            return true;
          }
        }

      } // end of case: q is neither horizontal nor vertical

    }

    // here and below both q and t are points

    if ( s.is_point() && r.is_point() && same_points(s, r) ) {
      // MK::ERROR: write this code using the compare_x_2 and
      //    compare_y_2 predicates instead of computing the inner
      //    product...
      // philaris: adaptation to Linf

      CGAL_SDG_DEBUG(std::cout <<
          "debug infcf s=r and is point" << std::endl;);

      Comparison_result cmpxst = cmpx(s.point(), t.point());
      Comparison_result cmpxtq = cmpx(t.point(), q.point());
      Comparison_result cmpyst = cmpy(s.point(), t.point());
      Comparison_result cmpytq = cmpy(t.point(), q.point());

      //Sign sgn1 = -sign_of( cmpxst * cmpxtq + cmpyst * cmpytq );
      Sign sgn1 = CGAL::compare(0, cmpxst * cmpxtq + cmpyst * cmpytq);

      CGAL_assertion( sgn1 != ZERO );

      CGAL_SDG_DEBUG(
          std::cout
          << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
          << q << ' ' << s << ' ' << r << ' ' << t
          << ' ' << sgn << " returns "
          << (sgn1 == POSITIVE) << std::endl;);
      return (sgn1 == POSITIVE);
    }

    // still here, both q and t are points

    if ( s.is_segment() && r.is_segment() && same_segments(s, r) ) {
      CGAL_SDG_DEBUG(std::cout <<
          "debug infcf s=r and is segment" << std::endl;);

      if ( same_points(q, s.source_site()) ||
           same_points(q, s.target_site())   )  {
        // here point q is endpoint of segment s
        Site_2 ss;
        if ( same_points(q, s.source_site()) ) {
          ss = s.target_site();
        } else {
          ss = s.source_site();
        }
        // MK::ERROR: write this code using the compare_x_2 and
        //    compare_y_2 predicates instead of computing the inner
        //    product...
        // philaris: adaptation to Linf

        Comparison_result cmpxst = cmpx(ss.point(), t.point());
        Comparison_result cmpxtq = cmpx(t.point(), q.point());
        Comparison_result cmpyst = cmpy(ss.point(), t.point());
        Comparison_result cmpytq = cmpy(t.point(), q.point());

        //Sign sgn1 = -sign_of( cmpxst * cmpxtq + cmpyst * cmpytq );
        Sign sgn1 = CGAL::compare(0, cmpxst * cmpxtq + cmpyst * cmpytq);

        CGAL_assertion( sgn1 != ZERO );

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << (sgn1 == POSITIVE) << std::endl;);

        return (sgn1 == POSITIVE);
      }
      else
      {
        // here point q is not endpoint of segment s

        CGAL_assertion(sgn == NEGATIVE);

        CGAL_assertion(! is_site_h_or_v(s));

        // compute infinite square with corner at q
        // and with center at infinity at
        // direction SE, NE, NW, or SW;
        // the direction goes from segment s to point q

        Line_2 l = compute_supporting_line(s.supporting_site());

        Point_2 phor =
          compute_horizontal_projection(l, q.point());

        Point_2 pver =
          compute_vertical_projection(l, q.point());

        // if (phor, q, pver) is a right turn,
        // then lines are (phor,q) and (q,pver)
        // else if (phor, q, pver) is a left turn,
        // then lines are (pver,q) and (q,phor)

        Line_2 l1, l2;

        if (CGAL::orientation(phor, q.point(), pver) == RIGHT_TURN) {
          l1 = compute_line_from_to(phor, q.point());
          l2 = compute_line_from_to(q.point(), pver);
        } else {
          l1 = compute_line_from_to(pver, q.point());
          l2 = compute_line_from_to(q.point(), phor);
        }

        // the square is defined by the intersection
        // of the positive halfspaces supported
        // by lines l1 and l2

        Oriented_side osl1 =
          oriented_side_of_line(l1, t.point());

        Oriented_side osl2 =
          oriented_side_of_line(l2, t.point());

        CGAL_SDG_DEBUG(std::cout << "debug iecf osl1=" << osl1
          << " osl2=" << osl2 << std::endl;);

        Boolean retval =
          (osl1 !=
           ON_NEGATIVE_SIDE) &&
          (osl2 !=
           ON_NEGATIVE_SIDE)     ;

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << retval << std::endl;);

        return retval;

      }
    }

    // philaris: here there is significant difference of Linf
    // with relation to L2

    // q is on the Linf convex hull with neighbors (in this hull)
    // s and r

    // here q is a point and s is different from r

    if (sgn == NEGATIVE) {
      CGAL_SDG_DEBUG(std::cout <<
          "debug infinite-edge-int-cf special NEG" << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug infcf special (q,s,r,t,sgn)= "
        << q << ' ' << s << ' ' << r << ' ' << t
        << ' ' << sgn << std::endl;);

      if (s.is_point() && r.is_point()) {
        if ((bounded_side_of_bbox(
                q.point(), s.point(), t.point()) == ON_BOUNDED_SIDE) &&
            (bounded_side_of_bbox(
                q.point(), r.point(), t.point()) == ON_BOUNDED_SIDE))
        {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << false << std::endl;);
          return false;
          // otherwise, it will return true later
        }
      } else {
        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << true << std::endl;);
        return true;
      }
    } else if (sgn == POSITIVE) {
      CGAL_SDG_DEBUG(
          std::cout << "debug infinite-edge-int-cf special POS"
          << std::endl;);
      if (s.is_point() && r.is_point()) {
        if ((bounded_side_of_bbox(
               t.point(), s.point(), q.point()) ==
                 ON_BOUNDED_SIDE) &&
            (bounded_side_of_bbox(
               t.point(), r.point(), q.point()) ==
                 ON_BOUNDED_SIDE))
        {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << true << std::endl;);
          return true;
          // otherwise it will return false later
        }
      } else if (s.is_segment() && r.is_segment()) {
        // here s and r are both segments

        // philaris: I assert that q is an endpoint of both s and r

        CGAL_assertion(
            same_points(q, s.source_site()) ||
            same_points(q, s.target_site())
            );

        CGAL_assertion(
            same_points(q, r.source_site()) ||
            same_points(q, r.target_site())
            );

        // Since s and r are different, the infinite vertices
        // (q,s,inf) and (q,inf,r) of the Voronoi diagram
        // have the following property:
        // if they are connected, their line has +-45 deg slope
        // and this line goes through point q.
        // This line is Linf-perpendicular to both s and r.
        // Since q is not in conflict with any of (q,s,inf),
        // (q,inf,r), it is enough to check if t is in the
        // infinite box with corner q containing neither
        // of these: (q,s,inf), (q,inf,r), any point of s, r.

        Point_2 otherpnt =
          (same_points(q, s.source_site())) ?
          s.segment().target():
          s.segment().source();

        CGAL_assertion(
            or_linf(t.point(), q.point(), otherpnt) ==
            DEGENERATE );

        Bounded_side bside =
          bounded_side_of_bbox(t.point(), otherpnt, q.point());

        CGAL_assertion(bside != ON_BOUNDARY);

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << (bside == ON_BOUNDED_SIDE) << std::endl;);
        return (bside == ON_BOUNDED_SIDE);

      } else {
        // here one of s, r is point and the other is a segment
        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn) "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " tocheck one seg one pnt in s, r"
            << std::endl; ) ;
      }
    }
    // here it might be sgn == ZERO

    CGAL_SDG_DEBUG(
        std::cout
        << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
        << q << ' ' << s << ' ' << r << ' ' << t
        << ' ' << sgn << " returns "
        << ( sgn == NEGATIVE ) << std::endl;);
    return ( sgn == NEGATIVE );
  }

};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H

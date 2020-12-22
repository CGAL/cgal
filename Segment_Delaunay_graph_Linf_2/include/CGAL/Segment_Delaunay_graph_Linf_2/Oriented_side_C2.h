// Copyright (c) 2015  Università della Svizzera italiana.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_x_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_y_2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//-----------------------------------------------------------------------------



template<class K, class Method_tag>
class Oriented_side_C2
  : public Basic_predicates_C2<K>
{
private:

  typedef Basic_predicates_C2<K>              Base;
  typedef Voronoi_vertex_C2<K,Method_tag>     Voronoi_vertex_2;

  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Line_2               Line_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;
  typedef typename Base::Comparison_result    Comparison_result;

  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K> Are_same_points_2;
  Are_same_points_2                same_points;

  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K> Are_same_segments_2;
  Are_same_segments_2              same_segments;

  typedef SegmentDelaunayGraph_2::Compare_x_2<K> Compare_x_2_Sites_Type;
  Compare_x_2_Sites_Type scmpx;
  typedef SegmentDelaunayGraph_2::Compare_y_2<K> Compare_y_2_Sites_Type;
  Compare_y_2_Sites_Type scmpy;

  using Base::compute_supporting_line;
  using Base::compute_linf_perpendicular;
  using Base::oriented_side_of_line;
  using Base::compute_vertical_projection;
  using Base::compute_horizontal_projection;
  using Base::is_site_h_or_v;
  using Base::compare_distance_to_point_linf;

public:
  typedef typename Base::Oriented_side        Oriented_side;
  typedef Oriented_side                       result_type;
  typedef Site_2                              argument_type;

  // computes the oriented side of the Voronoi vertex of s1, s2, s3
  // wrt the line that passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
                           const Site_2& s3,
                           const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    Voronoi_vertex_2 v(s1, s2, s3);
    Line_2 l = compute_supporting_line( s.supporting_site() );
    Line_2 lp = compute_linf_perpendicular(l, p.point());

    Oriented_side retval = v.oriented_side(lp);

    CGAL_SDG_DEBUG(std::cout
        << "debug: Oriented_side_C2 (s1,s2,s3,s,p)= ("
        << s1 << ") (" << s2 << ") (" << s3 << ") ("
        << s << ") (" << p << ") "
        << "returns " << retval
        << std::endl;);

    return retval;
  }

  // tie breaker for finite vertex
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
                           const Site_2& s3,
                           const Site_2& s, const Site_2& p,
                           const Point_2 & pt) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );
    CGAL_USE(pt);

    Voronoi_vertex_2 v(s1, s2, s3);
    Line_2 l = compute_supporting_line( s.supporting_site() );
    Line_2 lp = compute_linf_perpendicular(l, p.point());

    Oriented_side retval = v.oriented_side(lp);

    bool is_s1_pnt = s1.is_point();
    bool is_s2_pnt = s2.is_point();
    bool is_s3_pnt = s3.is_point();

    CGAL_assertion( ! (is_s1_pnt && is_s2_pnt && is_s3_pnt) );

    // tie breaker in case at least one site is a point and
    // the segment to split is not axis-parallel
    if ((is_s1_pnt || is_s2_pnt || is_s3_pnt) &&
        (! is_site_h_or_v(s))) {
      if (retval == ON_ORIENTED_BOUNDARY) {
        unsigned int num_pts =
          (is_s1_pnt? 1 : 0) +
          (is_s2_pnt? 1 : 0) +
          (is_s3_pnt? 1 : 0) ;
        CGAL_assertion( num_pts == 1 || num_pts == 2 );
        CGAL_SDG_DEBUG(std::cout
            << "debug: Oriented_side_C2 num_pts=" << num_pts
            << std::endl;);

        if (num_pts == 1) {
          const Site_2 * site_ptr;
          if (is_s1_pnt) {
            site_ptr = &s1;
          } else if (is_s2_pnt) {
            site_ptr = &s2;
          } else {
            site_ptr = &s3;
          }
          FT dist;
          bool is_cand = test_candidate(*site_ptr, p, v, dist);
          if (is_cand) {
            // override retval = 0
            retval = - oriented_side_of_line(lp, site_ptr->point());
          }
        } else { // num_pts == 2
          const Site_2 * a;
          const Site_2 * b;

          if (! is_s1_pnt) {
            a = &s2;
            b = &s3;
          } else if (! is_s2_pnt) {
            a = &s1;
            b = &s3;
          } else { // not is_s3_pnt
            a = &s1;
            b = &s2;
          }

          FT distpa (0);
          bool is_a_cand = test_candidate(*a, p, v, distpa);
          FT distpb (0);
          bool is_b_cand = test_candidate(*b, p, v, distpb);
          unsigned int num_candidates =
            (is_a_cand? 1 : 0) +
            (is_b_cand? 1 : 0) ;
          CGAL_SDG_DEBUG(std::cout
              << "debug: Oriented_side_C2 two points, num_candidates="
              << num_candidates << std::endl;);
          const Site_2 * test_site_ptr;
          if (num_candidates == 1) {
            test_site_ptr = (is_a_cand) ? a : b;
            retval = - oriented_side_of_line(lp, test_site_ptr->point());
          } else if (num_candidates == 2) {
            Comparison_result testab = CGAL::compare(distpa, distpb);
            if (testab != EQUAL) {
              test_site_ptr = (testab == SMALLER) ? a : b;
              retval = - oriented_side_of_line(lp, test_site_ptr->point());
            }
          } // end of case of num_candidates == 2
        } // end of case of num_pts == 2
      } // end of case ON_ORIENTED_BOUNDARY
    } // end of case of at least a point and non-hv seg s

    CGAL_SDG_DEBUG(std::cout
        << "debug: Oriented_side_C2 (s1,s2,s3,s,p)= ("
        << s1 << ") (" << s2 << ") (" << s3 << ") ("
        << s << ") (" << p << ") "
        << "tiebreaker returns " << retval
        << std::endl;);

    return retval;
  }

  // computes the oriented side of the Voronoi vertex of s1, s2, inf
  // wrt the line that passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
                           const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    Line_2 lseg = compute_supporting_line( s.supporting_site() );
    Line_2 lp = compute_linf_perpendicular(lseg, p.point());
    bool has_lseg_neg_slope;

    // Voronoi_vertex_2 v(s1, s2, inf);
    // compute linf projection of v(s1, s2, inf) on segment s,
    // which will be the test point for the orientation test

    Point_2 testpnt;

    bool is_s1_segment = s1.is_segment();
    bool is_s2_segment = s2.is_segment();

    CGAL_assertion(
        (is_s1_segment && (same_segments(s, s1) ||
                            same_segments(s, s1.supporting_site())))
        ||
        (is_s2_segment && (same_segments(s, s2) ||
                            same_segments(s, s2.supporting_site()))));

    bool are_both_segments = is_s1_segment && is_s2_segment;

    // boolean variable of:
    // point in {s1,s2} being endpoint of the segment in {s1,s2}
    bool are_endp_s1s2 =
           (is_s1_segment &&
            ( same_points(s2, s1.source_site()) ||
              same_points(s2, s1.target_site())   ) )
           ||
           (is_s2_segment &&
            ( same_points(s1, s2.source_site()) ||
              same_points(s1, s2.target_site())   ) )  ;

    if (are_both_segments) {
      // the two segments must have a common endpoint,
      // which is the linf projection

      CGAL_assertion( is_site_h_or_v(s1) || is_site_h_or_v(s2) );

      if (same_points(s1.source_site(), s2.source_site()) ||
          same_points(s1.source_site(), s2.target_site())   ) {
        testpnt = s1.source_site().point();
      } else {
        CGAL_assertion(
          same_points(s1.target_site(), s2.source_site()) ||
          same_points(s1.target_site(), s2.target_site())   );
        testpnt = s1.target_site().point();
      }

    } else {
      // here, there is a point and a segment in {s1, s2}

      if ( are_endp_s1s2 )
      { // here the point in {s1,s2}
        // is endpoint of the segment in {s1,s2}
        if (is_s1_segment) {
          testpnt = s2.point();
        } else {
          CGAL_assertion(is_s2_segment);
          testpnt = s1.point();
        }
      } // end of case: point is endpoint of segment
      else {
        // here, the point is not endpoint of the segment

        CGAL_SDG_DEBUG(std::cout << "debug: Oriented_side_C2 (s1,s2,s,p)= ("
              << s1 << ") (" << s2 << ") ("
              << s << ") (" << p << ") "
              << "case of s1/s2 no endpoint relation"
              << std::endl;);

        CGAL_assertion( ! is_site_h_or_v(s) );

        has_lseg_neg_slope =
          CGAL::sign(lseg.a()) == CGAL::sign(lseg.b());

        if (has_lseg_neg_slope) {
          if (is_s1_segment) {
            testpnt =
              compute_horizontal_projection(lseg, s2.point());
          } else {
            testpnt =
              compute_vertical_projection(lseg, s1.point());
          }
        } else {
          // here, segment has positive slope
          if (is_s1_segment) {
            testpnt =
              compute_vertical_projection(lseg, s2.point());
          } else {
            testpnt =
              compute_horizontal_projection(lseg, s1.point());
          }
        } // end of case: seg has positive slope
      } // end of case: point is not endpoint of segment

    } // end of case: a point and a segment in {s1, s2}


    Oriented_side retval =
      oriented_side_of_line(lp, testpnt);

    CGAL_SDG_DEBUG(std::cout << "debug: Oriented_side_C2 (s1,s2,s,p)= ("
              << s1 << ") (" << s2 << ") ("
              << s << ") (" << p << ") "
              << "returns " << retval
              << std::endl;);

    return retval;
  }

  // tie breaker for infinite vertex
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
                           const Site_2& s, const Site_2& p,
                           const Point_2 & pt) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );
    CGAL_USE(pt);

    Line_2 lseg = compute_supporting_line( s.supporting_site() );
    Line_2 lp = compute_linf_perpendicular(lseg, p.point());

    // Voronoi_vertex_2 v(s1, s2, inf);
    // compute linf projection of v(s1, s2, inf) on segment s,
    // which will be the test point for the orientation test

    Point_2 testpnt;

    bool is_s1_segment = s1.is_segment();
    bool is_s2_segment = s2.is_segment();

    CGAL_assertion(
        (is_s1_segment && (same_segments(s, s1) ||
                            same_segments(s, s1.supporting_site())))
        ||
        (is_s2_segment && (same_segments(s, s2) ||
                            same_segments(s, s2.supporting_site()))));

    bool are_both_segments = is_s1_segment && is_s2_segment;

    // boolean variable of:
    // point in {s1,s2} being endpoint of the segment in {s1,s2}
    bool are_endp_s1s2 =
           (is_s1_segment &&
            ( same_points(s2, s1.source_site()) ||
              same_points(s2, s1.target_site())   ) )
           ||
           (is_s2_segment &&
            ( same_points(s1, s2.source_site()) ||
              same_points(s1, s2.target_site())   ) )  ;

    if (are_both_segments) {
      // the two segments must have a common endpoint,
      // which is the linf projection

      CGAL_assertion( is_site_h_or_v(s1) || is_site_h_or_v(s2) );

      if (same_points(s1.source_site(), s2.source_site()) ||
          same_points(s1.source_site(), s2.target_site())   ) {
        testpnt = s1.source_site().point();
      } else {
        CGAL_assertion(
          same_points(s1.target_site(), s2.source_site()) ||
          same_points(s1.target_site(), s2.target_site())   );
        testpnt = s1.target_site().point();
      }

    } else {
      // here, there is a point and a segment in {s1, s2}

      if ( are_endp_s1s2 )
      { // here the point in {s1,s2}
        // is endpoint of the segment in {s1,s2}
        if (is_s1_segment) {
          testpnt = s2.point();
        } else {
          CGAL_assertion(is_s2_segment);
          testpnt = s1.point();
        }
      } // end of case: point is endpoint of segment
      else {
        // here, the point is not endpoint of the segment

        CGAL_SDG_DEBUG(std::cout << "debug: Oriented_side_C2 (s1,s2,s,p)= ("
              << s1 << ") (" << s2 << ") ("
              << s << ") (" << p << ") "
              << "case of s1/s2 no endpoint relation"
              << std::endl;);

        CGAL_assertion( ! is_site_h_or_v(s) );

        bool has_lseg_neg_slope =
          CGAL::sign(lseg.a()) == CGAL::sign(lseg.b());

        if (has_lseg_neg_slope) {
          if (is_s1_segment) {
            testpnt =
              compute_horizontal_projection(lseg, s2.point());
          } else {
            testpnt =
              compute_vertical_projection(lseg, s1.point());
          }
        } else {
          // here, segment has positive slope
          if (is_s1_segment) {
            testpnt =
              compute_vertical_projection(lseg, s2.point());
          } else {
            testpnt =
              compute_horizontal_projection(lseg, s1.point());
          }
        } // end of case: seg has positive slope
      } // end of case: point is not endpoint of segment

    } // end of case: a point and a segment in {s1, s2}


    Oriented_side retval =
      oriented_side_of_line(lp, testpnt);

    if (retval == ON_ORIENTED_BOUNDARY) {
      // philaris: tocheck this later
      CGAL_assertion(! are_both_segments);
      // philaris: tocheck this later
      CGAL_assertion(! are_endp_s1s2);

      CGAL_SDG_DEBUG(std::cout << "debug: Oriented_side_C2 (s1,s2,s,p)= ("
              << s1 << ") (" << s2 << ") ("
              << s << ") (" << p << ") "
              << "trying to fix ZERO"
              << std::endl;);
      if (is_s1_segment) {
        testpnt = s2.point();
      } else {
        testpnt = s1.point();
      }
      retval = - oriented_side_of_line(lp, testpnt);
    }


    CGAL_SDG_DEBUG(std::cout << "debug: Oriented_side_C2 (s1,s2,s,p)= ("
              << s1 << ") (" << s2 << ") ("
              << s << ") (" << p << ") "
              << "tiebreaker returns " << retval
              << std::endl;);

    return retval;
  }


private:

  inline
  bool test_candidate(const Site_2& c, const Site_2& p,
                      const Voronoi_vertex_2& v, FT& distpc)
    const
  {
    if (c.is_segment()) return false;
    // here c is point
    if (scmpx(p, c) == EQUAL) {
      distpc = CGAL::abs(p.point().y() - c.point().y());
      FT squareside = FT(2)*CGAL::abs(p.point().y() - v.point().y());
      if (CGAL::compare(distpc, squareside) == SMALLER) {
        return true;
      }
    } else if (scmpy(p, c) == EQUAL) {
      distpc = CGAL::abs(p.point().x() - c.point().x());
      FT squareside = FT(2)*CGAL::abs(p.point().x() - v.point().x());
      if (CGAL::compare(distpc, squareside) == SMALLER) {
        return true;
      }
    }
    return false;
  }

};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_C2_H

// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ORIENTED_SIDE_OF_BISECTOR_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_ORIENTED_SIDE_OF_BISECTOR_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//------------------------------------------------------------------------
//------------------------------------------------------------------------

template<class K, class Method_tag>
class Oriented_side_of_bisector_C2
  : public Basic_predicates_C2<K>
{
private:
  typedef Basic_predicates_C2<K>            Base;
  typedef typename Base::Line_2             Line_2;
  typedef typename Base::RT                 RT;
  typedef typename Base::FT                 FT;
  typedef typename Base::Comparison_result  Comparison_result;
  typedef typename Base::Sign               Sign;
  typedef typename Base::Orientation        Orientation;

  typedef Are_same_points_C2<K>             Are_same_points_2;
  typedef Are_same_segments_C2<K>           Are_same_segments_2;

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

public:
  typedef typename Base::Oriented_side      Oriented_side;
  typedef typename Base::Site_2             Site_2;
  typedef typename Base::Point_2            Point_2;
  typedef typename Base::Segment_2          Segment_2;

  using Base::compute_supporting_line;
  using Base::compute_perpendicular;
  using Base::opposite_line;
  using Base::oriented_side_of_line;
  using Base::compute_squared_distance;
private:
  bool is_endpoint(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    return
      same_points(p, s.source_site()) || same_points(p, s.target_site());
  }

  bool is_degenerate(const Site_2& s) const
  {
    CGAL_precondition( s.is_segment() );
    return same_points( s.source_site(), s.target_site() );
  }

  //-----------------------------------------------------------------

  Comparison_result
  compare_distances_pp(const Site_2& p1, const Site_2& p2,
                       const Site_2& q) const
  {
    CGAL_precondition( p1.is_point() && p2.is_point() );
    CGAL_precondition( !same_points(p1,p2) );

    if ( same_points(q, p1) ) { return SMALLER; }
    if ( same_points(q, p2) ) { return LARGER; }

    return
      compare_distance_to_point(q.point(), p1.point(), p2.point());
  }

  //-----------------------------------------------------------------

  Comparison_result
  compare_distances_sp(const Site_2& s, const Site_2& p,
                       const Site_2& q) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );
    CGAL_precondition( !is_degenerate(s) );

    if ( same_points(q, p) ) { return LARGER; }
    if ( same_points(q, s.source_site()) ) { return SMALLER; }
    if ( same_points(q, s.target_site()) ) { return SMALLER; }


    bool is_src = same_points(p, s.source_site());
    bool is_trg = same_points(p, s.target_site());

    if ( is_src || is_trg ) {
      Line_2 ls = compute_supporting_line(s.supporting_site());
      Line_2 lp = compute_perpendicular(ls, p.point());

      if ( is_trg ) {
        lp = opposite_line(lp);
      }

      Oriented_side os = oriented_side_of_line(lp, q.point());

      if ( os == ON_POSITIVE_SIDE ) {
        return LARGER;
      } else if ( os == ON_NEGATIVE_SIDE) {
        return SMALLER;
      } else {
        return EQUAL;
      }
    }

    Point_2 pp = p.point(), qq = q.point();

    RT d2_p = compute_squared_distance(pp, qq);

    Point_2 ssrc = s.source(), strg = s.target();

    Line_2 ls = compute_supporting_line(s.supporting_site());
    Line_2 lsrc = compute_perpendicular(ls, ssrc);
    Line_2 ltrg = compute_perpendicular(ls, strg);

    Oriented_side os_src = oriented_side_of_line(lsrc, qq);
    if ( os_src != ON_NEGATIVE_SIDE ) {
      RT d2_s = compute_squared_distance(qq, ssrc);
      return CGAL::compare(d2_s, d2_p);
    }

    Oriented_side os_trg = oriented_side_of_line(ltrg, qq);
    if ( os_trg != ON_POSITIVE_SIDE ) {
      RT d2_s = compute_squared_distance(qq, strg);
      return CGAL::compare(d2_s, d2_p);
    }

    std::pair<RT,RT> d2_s = compute_squared_distance(qq, ls);
    return CGAL::compare(d2_s.first, d2_p * d2_s.second);
  }

  //-----------------------------------------------------------------

  Comparison_result
  compare_distances_ss(const Site_2& s1, const Site_2& s2,
                       const Site_2& q) const
  {
    CGAL_precondition( s1.is_segment() && s2.is_segment() );
    CGAL_precondition( !is_degenerate(s1) );
    CGAL_precondition( !is_degenerate(s2) );

    bool is_on_s1 = is_endpoint(q, s1);
    bool is_on_s2 = is_endpoint(q, s2);

    if ( is_on_s1 && is_on_s2 ) {
      return EQUAL;
    } else if ( is_on_s1 && !is_on_s2 ) {
      return SMALLER;
    } else if ( !is_on_s1 && is_on_s2 ) {
      return LARGER;
    }

    if ( same_segments(s1, s2) ) {
      return EQUAL;
    }


    Point_2 qq = q.point();

    Point_2 ssrc1 = s1.source(), strg1 = s1.target();

    Line_2 ls1 = compute_supporting_line(s1.supporting_site());
    Line_2 lsrc1 = compute_perpendicular(ls1, ssrc1);
    Line_2 ltrg1 = compute_perpendicular(ls1, strg1);

    Point_2 ssrc2 = s2.source(), strg2 = s2.target();

    Line_2 ls2 = compute_supporting_line(s2.supporting_site());
    Line_2 lsrc2 = compute_perpendicular(ls2, ssrc2);
    Line_2 ltrg2 = compute_perpendicular(ls2, strg2);

    // idx1 and idx2 indicate if q is to the left (source endpoint
    // side), the right side (target endpoint side) or inside
    // the band of s1 and s2 respectively; if q is on the boundary of
    // the band we assign it to the adjacent halfplane; for left
    // halfplane the value is -1; for the band the value is 0; for the
    // right halfplane the value is 1.
    int idx1(0), idx2(0);

    Oriented_side os_src1 = oriented_side_of_line(lsrc1, qq);
    if ( os_src1 != ON_NEGATIVE_SIDE ) {
      idx1 = -1;
    } else {
      Oriented_side os_trg1 = oriented_side_of_line(ltrg1, qq);
      if ( os_trg1 != ON_POSITIVE_SIDE ) {
        idx1 = 1;
      }
    }

    Oriented_side os_src2 = oriented_side_of_line(lsrc2, qq);
    if ( os_src2 != ON_NEGATIVE_SIDE ) {
      idx2 = -1;
    } else {
      Oriented_side os_trg2 = oriented_side_of_line(ltrg2, qq);
      if ( os_trg2 != ON_POSITIVE_SIDE ) {
        idx2 = 1;
      }
    }

    CGAL_assertion( idx1 >= -1 && idx1 <= 1 );
    CGAL_assertion( idx2 >= -1 && idx2 <= 1 );

    if ( idx1 == -1 ) {
      RT d2_s1 = compute_squared_distance(qq, ssrc1);
      if ( idx2 == -1 ) {
        if ( same_points(s1.source_site(), s2.source_site()) ) {
          return EQUAL;
        }
        RT d2_s2 = compute_squared_distance(qq, ssrc2);
        return CGAL::compare(d2_s1, d2_s2);
      } else if ( idx2 == 1 ) {
        if ( same_points(s1.source_site(), s2.target_site()) ) {
          return EQUAL;
        }
        RT d2_s2 = compute_squared_distance(qq, strg2);
        return CGAL::compare(d2_s1, d2_s2);
      } else {
        std::pair<RT,RT> d2_s2 = compute_squared_distance(qq, ls2);
        return CGAL::compare(d2_s1 * d2_s2.second, d2_s2.first);
      }
    } else if ( idx1 == 1 ) {
      RT d2_s1 = compute_squared_distance(qq, strg1);
      if ( idx2 == -1 ) {
        if ( same_points(s1.target_site(), s2.source_site()) ) {
          return EQUAL;
        }
        RT d2_s2 = compute_squared_distance(qq, ssrc2);
        return CGAL::compare(d2_s1, d2_s2);
      } else if ( idx2 == 1 ) {
        if ( same_points(s1.target_site(), s2.target_site()) ) {
          return EQUAL;
        }
        RT d2_s2 = compute_squared_distance(qq, strg2);
        return CGAL::compare(d2_s1, d2_s2);
      } else {
        std::pair<RT,RT> d2_s2 = compute_squared_distance(qq, ls2);
        return CGAL::compare(d2_s1 * d2_s2.second, d2_s2.first);
      }
    }

    CGAL_assertion( idx1 == 0 );
    std::pair<RT,RT> d2_s1 = compute_squared_distance(qq, ls1);
    if ( idx2 == -1 ) {
      RT d2_s2 = compute_squared_distance(qq, ssrc2);
      return CGAL::compare(d2_s1.first, d2_s2 * d2_s1.second);
    } else if ( idx2 == 1 ) {
      RT d2_s2 = compute_squared_distance(qq, strg2);
      return CGAL::compare(d2_s1.first, d2_s2 * d2_s1.second);
    }

    CGAL_assertion( idx2 == 0 );
    std::pair<RT,RT> d2_s2 = compute_squared_distance(qq, ls2);
    return CGAL::compare(d2_s1.first * d2_s2.second,
                         d2_s2.first * d2_s1.second);
  }

  //-----------------------------------------------------------------
  //-----------------------------------------------------------------

  Oriented_side
  compare_distances(const Site_2& t1, const Site_2& t2,
                    const Site_2& q) const
  {
    Comparison_result r;

    if ( t1.is_point() && t2.is_point() ) {
      r = compare_distances_pp(t1, t2, q);
    } else if ( t1.is_segment() && t2.is_point() ) {
      r = compare_distances_sp(t1, t2, q);
    } else if ( t1.is_point() && t2.is_segment() ) {
      r = opposite( compare_distances_sp(t2, t1, q) );
    } else {
      r = compare_distances_ss(t1, t2, q);
    }

    return -r;
  }

public:
  typedef Oriented_side              result_type;
  typedef Site_2                     argument_type;


  Oriented_side
  operator()(const Site_2& t1, const Site_2& t2, const Site_2& q) const
  {
    CGAL_precondition( q.is_point() );
    return compare_distances(t1, t2, q);
  }
};

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ORIENTED_SIDE_OF_BISECTOR_C2_H

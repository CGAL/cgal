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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_OF_BISECTOR_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_OF_BISECTOR_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//------------------------------------------------------------------------
//------------------------------------------------------------------------

template<class K, class Method_tag>
class Oriented_side_of_bisector_C2
  : public Basic_predicates_C2<K>
{
private:
  typedef Basic_predicates_C2<K>             Base;
  typedef typename Base::Line_2              Line_2;
  typedef typename Base::RT                  RT;
  typedef typename Base::FT                  FT;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Sign                Sign;
  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>   Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K> Are_same_segments_2;

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

public:
  typedef typename Base::Oriented_side      Oriented_side;
  typedef typename Base::Site_2             Site_2;
  typedef typename Base::Point_2            Point_2;
  typedef typename Base::Segment_2          Segment_2;

  using Base::compute_supporting_line;
  using Base::compute_linf_perpendicular;
  using Base::compute_linf_projection_hom;
  using Base::opposite_line;
  using Base::oriented_side_of_line;
  using Base::compute_linf_distance;
  using Base::compare_distance_to_point_linf;

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

    //CGAL_SDG_DEBUG(std::cout << "debug compare_distances_pp" << std::endl;);

    return
      compare_distance_to_point_linf(q.point(), p1.point(), p2.point());
  }

  //-----------------------------------------------------------------

  Comparison_result
  compare_distances_sp(const Site_2& s, const Site_2& p,
		       const Site_2& q) const
  {

    //CGAL_SDG_DEBUG(std::cout << "debug compare_distances_sp entering "
    //  << "(s =" << s << ") p=(" << p << ") q=(" << q << ")"
    //  << std::endl;);

    CGAL_precondition( s.is_segment() && p.is_point() );
    CGAL_precondition( !is_degenerate(s) );

    if ( same_points(q, p) ) { return LARGER; }
    if ( same_points(q, s.source_site()) ) { return SMALLER; }
    if ( same_points(q, s.target_site()) ) { return SMALLER; }


    const bool is_src = same_points(p, s.source_site());
    const bool is_trg = same_points(p, s.target_site());

    if ( is_src || is_trg ) {
      Line_2 ls = compute_supporting_line(s.supporting_site());
      Line_2 lp = compute_linf_perpendicular(ls, p.point());

      if ( is_trg ) {
	lp = opposite_line(lp);
      }

      const Oriented_side os = oriented_side_of_line(lp, q.point());

      CGAL_SDG_DEBUG(std::cout << "debug compare_distances_sp "
          << " is_src=" << is_src << " is_trg=" << is_trg
          << " has os=" << os << std::endl;);

      if ( os == ON_POSITIVE_SIDE ) {
	return LARGER;
      } else if ( os == ON_NEGATIVE_SIDE) {
	return SMALLER;
      } else {
        // philaris: here, point is closer, not interior of segment
        //return LARGER;
        return EQUAL;
      }
    }

    const Point_2 pp = p.point(), qq = q.point();
    const Line_2 ls = compute_supporting_line(s.supporting_site());

    const Point_2 ssrc = s.source();
    const Line_2 lsrc = compute_linf_perpendicular(ls, ssrc);
    const Oriented_side os_src = oriented_side_of_line(lsrc, qq);
    if ( os_src != ON_NEGATIVE_SIDE ) {
      return compare_distance_to_point_linf(qq, ssrc, pp);
    }

    const Point_2 strg = s.target();
    const Line_2 ltrg = compute_linf_perpendicular(ls, strg);
    const Oriented_side os_trg = oriented_side_of_line(ltrg, qq);
    if ( os_trg != ON_POSITIVE_SIDE ) {
      return compare_distance_to_point_linf(qq, strg, pp);
    }

    const RT d2_p = compute_linf_distance(pp, qq);
    const std::pair<RT,RT> d2_s = compute_linf_distance(qq, ls);
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

    CGAL_SDG_DEBUG(std::cout << "debug compare_distances_ss "
        << "entering with s1=" << s1 << " s2=" << s2
        << " q=" << q << std::endl;);

    const bool is_on_s1 = is_endpoint(q, s1);
    const bool is_on_s2 = is_endpoint(q, s2);

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

    const Point_2 qq = q.point();

    const Point_2 ssrc1 = s1.source(), strg1 = s1.target();

    const Line_2 ls1 = compute_supporting_line(s1.supporting_site());
    const Line_2 lsrc1 = compute_linf_perpendicular(ls1, ssrc1);
    const Line_2 ltrg1 = compute_linf_perpendicular(ls1, strg1);

    const Point_2 ssrc2 = s2.source(), strg2 = s2.target();

    const Line_2 ls2 = compute_supporting_line(s2.supporting_site());
    const Line_2 lsrc2 = compute_linf_perpendicular(ls2, ssrc2);
    const Line_2 ltrg2 = compute_linf_perpendicular(ls2, strg2);

    // idx1 and idx2 indicate if q is to the left (source endpoint
    // side), the right side (target endpoint side) or inside
    // the band of s1 and s2 respectively; if q is on the boundary of
    // the band we assign it to the adjacent halfplane; for left
    // halfplane the value is -1; for the band the value is 0; for the
    // right halfplane the value is 1.
    int idx1(0), idx2(0);

    const Oriented_side os_src1 = oriented_side_of_line(lsrc1, qq);
    Oriented_side os_trg1;
    if ( os_src1 != ON_NEGATIVE_SIDE ) {
      idx1 = -1;
      os_trg1 = ON_POSITIVE_SIDE;
    } else {
      os_trg1 = oriented_side_of_line(ltrg1, qq);
      if ( os_trg1 != ON_POSITIVE_SIDE ) {
	idx1 = 1;
      }
    }

    const Oriented_side os_src2 = oriented_side_of_line(lsrc2, qq);
    Oriented_side os_trg2;
    if ( os_src2 != ON_NEGATIVE_SIDE ) {
      idx2 = -1;
      os_trg2 = ON_POSITIVE_SIDE;
    } else {
      os_trg2 = oriented_side_of_line(ltrg2, qq);
      if ( os_trg2 != ON_POSITIVE_SIDE ) {
	idx2 = 1;
      }
    }

    CGAL_SDG_DEBUG(std::cout << "debug compare_distances_ss "
        << " os_src1=" << os_src1 << " os_trg1=" << os_trg1
        << " os_src2=" << os_src2 << " os_trg2=" << os_trg2 << std::endl;);

    CGAL_SDG_DEBUG(std::cout << "debug compare_distances_ss "
        << " idx=" << idx1 << " idx2=" << idx2 << std::endl;);

    CGAL_assertion( idx1 >= -1 && idx1 <= 1 );
    CGAL_assertion( idx2 >= -1 && idx2 <= 1 );

    if ( idx1 == -1 ) {
      if ( idx2 == -1 ) {
	if ( same_points(s1.source_site(), s2.source_site()) ) {
	  return EQUAL;
	}
	return compare_distance_to_point_linf(qq, ssrc1, ssrc2);
      } else if ( idx2 == 1 ) {
	if ( same_points(s1.source_site(), s2.target_site()) ) {
	  return EQUAL;
	}
	return compare_distance_to_point_linf(qq, ssrc1, strg2);
      } else {
        const RT d2_s1 = compute_linf_distance(ssrc1, qq);
        const std::pair<RT,RT> d2_s2 = compute_linf_distance(qq, ls2);
        return CGAL::compare(d2_s1 * d2_s2.second, d2_s2.first);
      }
    } else if ( idx1 == 1 ) {
      if ( idx2 == -1 ) {
	if ( same_points(s1.target_site(), s2.source_site()) ) {
	  return EQUAL;
	}
	return compare_distance_to_point_linf(qq, strg1, ssrc2);
      } else if ( idx2 == 1 ) {
	if ( same_points(s1.target_site(), s2.target_site()) ) {
	  return EQUAL;
	}
	return compare_distance_to_point_linf(qq, strg1, strg2);
      } else {
        // here closest2 is inside segment s2
        const RT d2_s1 = compute_linf_distance(strg1, qq);
        const std::pair<RT,RT> d2_s2 = compute_linf_distance(qq, ls2);
        return CGAL::compare(d2_s1 * d2_s2.second, d2_s2.first);
      }
    } else {
      // here closest1 is inside segment s1
      CGAL_assertion( idx1 == 0 );
      const std::pair<RT,RT> d2_s1 = compute_linf_distance(qq, ls1);
      if ( idx2 == -1 ) {
        if (is_endpoint(s2.source_site(), s1)) { return SMALLER; }
        const RT d2_s2 = compute_linf_distance(qq, ssrc2);
        return CGAL::compare(d2_s1.first, d2_s2 * d2_s1.second);
      } else if ( idx2 == 1 ) {
        if (is_endpoint(s2.target_site(), s1)) { return SMALLER; }
        const RT d2_s2 = compute_linf_distance(qq, strg2);
        return CGAL::compare(d2_s1.first, d2_s2 * d2_s1.second);
      } else {
        CGAL_assertion( idx2 == 0 );
        const std::pair<RT,RT> d2_s2 = compute_linf_distance(qq, ls2);
        return CGAL::compare(d2_s1.first * d2_s2.second,
                             d2_s2.first * d2_s1.second);
      }
    }
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
      //CGAL_SDG_DEBUG(std::cout << "debug compare_distances pp result"
      //          << " t1=" << t1 << " t2=" << t2
      //          << " q=" << q << "  res=" << r << std::endl;);
    } else if ( t1.is_segment() && t2.is_point() ) {
      r = compare_distances_sp(t1, t2, q);
      //CGAL_SDG_DEBUG(std::cout << "debug compare_distances sp result"
      //          << " t1=" << t1 << " t2=" << t2
      //          << " q=" << q << "  res=" << r << std::endl;);
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

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_OF_BISECTOR_C2_H

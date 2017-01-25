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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_FINITE_EDGE_INTERIOR_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_FINITE_EDGE_INTERIOR_CONFLICT_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4800)
// complaint about performance where we can't do anything
#endif

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//------------------------------------------------------------------------

template<class K, class Method_tag>
class Finite_edge_interior_conflict_C2
  : public Basic_predicates_C2<K>
{
public:

  typedef Basic_predicates_C2<K>              Base;
  typedef Voronoi_vertex_C2<K,Method_tag>     Voronoi_vertex_2;

  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Line_2               Line_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;

  typedef typename Base::Comparison_result    Comparison_result;
  typedef typename Base::Sign                 Sign;
  typedef typename Base::Orientation          Orientation;
  typedef typename Base::Oriented_side        Oriented_side;
  typedef typename Base::Boolean              Boolean;

  typedef typename Base::Homogeneous_point_2  Homogeneous_point_2;


  using Base::opposite_line;
  using Base::compute_supporting_line;
  using Base::oriented_side_of_line;
  using Base::compare_linf_distances_to_line;
  using Base::compare_linf_distances_to_lines;
  using Base::compute_linf_perpendicular;
  //using Base::projection_on_line;
  using Base::compute_linf_projection_hom;
  using Base::compute_linf_projection_nonhom;
  using Base::compute_vertical_projection;
  using Base::compute_horizontal_projection;
  using Base::intersects_segment_interior_bbox;
  using Base::intersects_segment_positive_of_wedge;
  using Base::intersects_segment_negative_of_wedge;
  using Base::intersects_segment_interior_inf_box;
  using Base::intersects_segment_interior_inf_wedge_sp;
  using Base::is_site_h_or_v;
  using Base::zero_voronoi_area;


private:
  typedef typename Base::Compare_x_2          Compare_x_2_points;
  typedef typename Base::Compare_y_2          Compare_y_2_points;

  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>   Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K>
          Are_same_segments_2;

  typedef typename K::Intersections_tag       ITag;

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

  Compare_x_2_points cmpx;
  Compare_y_2_points cmpy;

private:

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------

  Boolean
  is_interior_in_conflict_both(const Site_2& p, const Site_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag tag) const
  {
    Boolean   in_conflict(false);

    if ( zero_voronoi_area(p, r, s) || zero_voronoi_area(q, r, s) ) {
      return true;
    }

    if ( p.is_point() && q.is_point() ) {
      in_conflict = is_interior_in_conflict_both_pp(p, q, r, s, t, tag);

    } else if ( p.is_segment() && q.is_segment() ) {

      in_conflict = is_interior_in_conflict_both_ss(p, q, r, s, t, tag);

    } else if ( p.is_point() && q.is_segment() ) {

      in_conflict = is_interior_in_conflict_both_ps(p, q, r, s, t, tag);

    } else { // p is a segment and q is a point

      in_conflict = is_interior_in_conflict_both_sp(p, q, r, s, t, tag);
    }

    return in_conflict;
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both_pp(const Site_2& sp, const Site_2& sq,
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag ) const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() );

    Point_2 p = sp.point(), q = sq.point();

    if ( t.is_point() ) { return true; }

    // here t is a segment
    Line_2 lt = compute_supporting_line(t.supporting_site());

    Oriented_side op, oq;

    if ( same_points(sp, t.source_site()) ||
	 same_points(sp, t.target_site()) ) {
      op = ON_ORIENTED_BOUNDARY;
    } else {
      op = oriented_side_of_line(lt, p);
    }

    if ( same_points(sq, t.source_site()) ||
	 same_points(sq, t.target_site()) ) {
      oq = ON_ORIENTED_BOUNDARY;
    } else {
      oq = oriented_side_of_line(lt, q);
    }


    if ((op == ON_POSITIVE_SIDE && oq == ON_NEGATIVE_SIDE) ||
	(op == ON_NEGATIVE_SIDE && oq == ON_POSITIVE_SIDE) ||
	(op == ON_ORIENTED_BOUNDARY || oq == ON_ORIENTED_BOUNDARY)) {
      return true;
    }

    Comparison_result res =
      compare_linf_distances_to_line(lt, p, q);

    if ( res == EQUAL ) { return true; }

    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);


    Line_2 lperp;
    if ( res == SMALLER ) {
      // p is closer to lt than q
      lperp = compute_linf_perpendicular(lt, p);
    } else {
      // q is closer to lt than p
      lperp = compute_linf_perpendicular(lt, q);
    }

    Oriented_side opqr = vpqr.oriented_side(lperp);
    Oriented_side oqps = vqps.oriented_side(lperp);

    return ( opqr == oqps );
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both_ss(const Site_2& p, const Site_2& q,
				  const Site_2& , const Site_2& ,
				  const Site_2& , Method_tag) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() );
    CGAL_USE(p);
    CGAL_USE(q);
    return true;
  }

  //--------------------------------------------------------------------

  Boolean
  is_interior_in_conflict_both_ps(const Site_2& p, const Site_2& q,
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag tag) const
  {
    CGAL_precondition( p.is_point() && q.is_segment() );

    if ( same_points(p, q.source_site()) ||
	 same_points(p, q.target_site()) ) {
      return false;
    }

    if ( t.is_point() ) {
      return is_interior_in_conflict_both_ps_p(p, q, r, s, t, tag);
    }
    return is_interior_in_conflict_both_ps_s(p, q, r, s, t, tag);
  }

  //--------------------------------------------------------------------

  Boolean
  is_interior_in_conflict_both_ps_p(const Site_2& p, const Site_2& q,
				    const Site_2& r, const Site_2& s,
				    const Site_2& t, Method_tag ) const
  {
    CGAL_precondition( t.is_point() );

    Line_2 lq = compute_supporting_line(q.supporting_site());

    Comparison_result res =
      compare_linf_distances_to_line(lq, p.point(), t.point());

    //if ( res != SMALLER ) { return true; }
    if (certainly( res != SMALLER ) ) { return true; }
    if (! is_certain( res != SMALLER ) ) {
      return indeterminate<Boolean>();
    }

    Voronoi_vertex_2 vpqr(p, q, r);
    Voronoi_vertex_2 vqps(q, p, s);

    Line_2 lperp = compute_linf_perpendicular(lq, p.point());

    Oriented_side opqr = vpqr.oriented_side(lperp);
    Oriented_side oqps = vqps.oriented_side(lperp);

    return (opqr == oqps);
  }

  //--------------------------------------------------------------------

  bool check_if_exact(const Site_2&, const Tag_false&) const
  {
    return true;
  }

  bool check_if_exact(const Site_2& t1, const Tag_true&) const
  {
    return t1.is_input();
  }

  Boolean
  is_interior_in_conflict_both_ps_s(const Site_2& sp, const Site_2& sq,
				    const Site_2& r, const Site_2& s,
				    const Site_2& st, Method_tag ) const
  {
    CGAL_precondition( st.is_segment() );

    Point_2 p = sp.point();
    //    Segment_2 q = sq.segment(), t = st.segment();

    Line_2 lq = compute_supporting_line(sq.supporting_site());

    if ( oriented_side_of_line(lq, p) == ON_NEGATIVE_SIDE ) {
      lq = opposite_line(lq);
    }

    if ( same_points(sp, st.source_site()) ||
	 same_points(sp, st.target_site()) ) {

      Line_2 lqperp = compute_linf_perpendicular(lq, p);

      Voronoi_vertex_2 vpqr(sp, sq, r);
      Voronoi_vertex_2 vqps(sq, sp, s);

      Oriented_side opqr = vpqr.oriented_side(lqperp);
      Oriented_side oqps = vqps.oriented_side(lqperp);

      CGAL_SDG_DEBUG( std::cout
          << "debug endpt case computing ondifparabarcs" << std::endl; );
      Boolean   on_different_parabola_arcs =
	 ((opqr == ON_NEGATIVE_SIDE) & (oqps == ON_POSITIVE_SIDE)) |
	 ((opqr == ON_POSITIVE_SIDE) & (oqps == ON_NEGATIVE_SIDE));

      //if ( !on_different_parabola_arcs ) { return true; }
      if (certainly( !on_different_parabola_arcs ) ) { return true; }
      if (! is_certain( !on_different_parabola_arcs ) ) {
        return indeterminate<Boolean>();
      }

      Site_2 t1;
      if ( same_points(sp, st.source_site()) ) {
	t1 = st.target_site();
      } else {
	t1 = st.source_site();
      }

      Oriented_side o_t1;

      if ( same_points(t1, sq.source_site()) ||
	   same_points(t1, sq.target_site()) ) {
	o_t1 = ON_ORIENTED_BOUNDARY;
      } else if (  !check_if_exact(t1, ITag()) &&
		   ( same_segments(t1.supporting_site(0),
				   sq.supporting_site()) ||
		     same_segments(t1.supporting_site(1),
				   sq.supporting_site()) )  ) {
	o_t1 = ON_ORIENTED_BOUNDARY;
      } else {
	o_t1 = oriented_side_of_line(lq, t1.point());
      }

      if ( o_t1 == ON_NEGATIVE_SIDE ) {
	return true;
      }

      Comparison_result res =
	compare_linf_distances_to_line(lq, p, t1.point());

      return ( res == LARGER );
    } // end of case: point p is an endpoint of segment t

    Line_2 lt = compute_supporting_line(st.supporting_site());

    if ( oriented_side_of_line(lt, p) == ON_NEGATIVE_SIDE ) {
      lt = opposite_line(lt);
    }

    Comparison_result res =
      CGAL::compare(lt.a() * lq.b(), lt.b() * lq.a());
    bool are_parallel = (res == EQUAL);

    if ( are_parallel ) {
      Sign sgn = CGAL::sign(lt.a() * lq.a() + lt.b() * lq.b());
      bool have_opposite_directions = (sgn == NEGATIVE);
      if ( have_opposite_directions ) { lq = opposite_line(lq); }

      if ( oriented_side_of_line(lq, p) == oriented_side_of_line(lt, p) ) {
	return true;
      }

      if ( have_opposite_directions ) {
	lq = opposite_line(lq);
      }
    }

    Line_2 l = compute_linf_perpendicular(lt, p);

    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    Oriented_side o_l_pqr = vpqr.oriented_side(l);
    Oriented_side o_l_qps = vqps.oriented_side(l);
    if (certainly( (o_l_pqr == ON_POSITIVE_SIDE) &
	           (o_l_qps == ON_NEGATIVE_SIDE) ) )
        return false;
    if (certainly( (o_l_pqr == ON_NEGATIVE_SIDE) &
	           (o_l_qps == ON_POSITIVE_SIDE) ) )
	return true;
    if (! is_certain((o_l_pqr == -o_l_qps) & (o_l_pqr != ZERO)))
        return indeterminate<Boolean>();

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>> HERE I NEED TO CHECK THE BOUNDARY CASES <<<<<<
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line_2 lqperp = compute_linf_perpendicular(lq, p);

    Oriented_side opqr = vpqr.oriented_side(lqperp);
    Oriented_side oqps = vqps.oriented_side(lqperp);

    CGAL_SDG_DEBUG( std::cout
        << "debug computing ondifparabarcs" << std::endl; );

    Boolean   on_different_parabola_arcs = (opqr == -oqps) & (opqr != ZERO);

    // if ( !on_different_parabola_arcs ) { return true; }
    if (certainly( !on_different_parabola_arcs ) ) { return true; }
    if (! is_certain( !on_different_parabola_arcs ) ) {
      return indeterminate<Boolean>(); }

    Homogeneous_point_2 pv = compute_linf_projection_hom(lq, p);
    Homogeneous_point_2 hp(p);

    pv = Base::midpoint(pv, hp);

    Oriented_side o_l_pv = oriented_side_of_line(l, pv);

    CGAL_assertion( o_l_pv != ON_ORIENTED_BOUNDARY );

    CGAL_assertion( o_l_pqr != ON_ORIENTED_BOUNDARY ||
		    o_l_qps != ON_ORIENTED_BOUNDARY );

    if ( o_l_pqr == ON_ORIENTED_BOUNDARY ) {
      return ( o_l_qps == o_l_pv );
    } else {
      return ( o_l_pqr == o_l_pv );
    }

  }

  //--------------------------------------------------------------------

  Boolean
  is_interior_in_conflict_both_sp(const Site_2& p, const Site_2& q,
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag tag) const
  {
    return is_interior_in_conflict_both_ps(q, p, s, r, t, tag);
  }


  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_touch(const Site_2& p, const Site_2& q,
				const Site_2& r, const Site_2& s,
				const Site_2& t, Method_tag tag) const
  {
    // checks if interior of voronoi edge is in conflict if both extrema
    // of the voronoi edge touch the corresponding circles.
    // return true if interior is in conflict; false otherwise
    if ( t.is_segment() ) { return false; }

#if 1
    CGAL_assertion( p.is_segment() || q.is_segment() );

    bool are_same_sites_rs (false);
    if (r.is_point() && s.is_point()) {
      if (same_points(r, s)) {
        are_same_sites_rs = true;
      }
    } else {
      if (r.is_segment() && s.is_segment()) {
        if (same_segments(r, s)) {
          are_same_sites_rs = true;
        }
      }
    }

    if (are_same_sites_rs) {
      // philaris: the text would give false
      // philaris: for segments, we have to refine later
      if (r.is_point() && s.is_point()) {
        CGAL_SDG_DEBUG(std::cout << "debug is_interior_in_conflict_touch"
            << " about to return false with same points" << std::endl; );
        return false;
      }
    } else {
      Voronoi_vertex_2 vpqr(p, q, r);
      Voronoi_vertex_2 vqps(q, p, s);

      if ( vpqr.incircle_no_easy(s) == ZERO &&
          vqps.incircle_no_easy(r) == ZERO ) {
        CGAL_SDG_DEBUG(std::cout << "debug is_interior_in_conflict_touch"
            << " about to return false with both ZERO" << std::endl; );
        return false;
      }
    }

    if ( p.is_segment() && q.is_segment() ) {
      return true;
    }
#else
    // OLD CODE: buggy if the edge is degenerate
    if ( (p.is_point() && q.is_point()) ||
	 (p.is_segment() && q.is_segment()) ) {
      return true;
    }
#endif

    if ( p.is_point() && q.is_segment() ) {
      Line_2 lq = compute_supporting_line(q.supporting_site());

      Comparison_result res =
	compare_linf_distances_to_line(lq, p.point(), t.point());

      return (res != SMALLER);
    }

    CGAL_SDG_DEBUG(std::cout << "debug is_interior_in_conflict_touch"
        << " recursion with permuted input" << std::endl; );
    return is_interior_in_conflict_touch(q, p, s, r, t, tag);
  }

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------


  bool
  is_interior_in_conflict_none(const Site_2& p, const Site_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag tag) const
  {
    if ( t.is_segment() ) { return false; }

    bool in_conflict(false);

    if ( p.is_point() && q.is_point() ) {
      in_conflict = is_interior_in_conflict_none_pp(p, q, r, s, t, tag);
    } else if ( p.is_point() && q.is_segment() ) {
      in_conflict = is_interior_in_conflict_none_ps(p, q, r, s, t, tag);
    } else if ( p.is_segment() && q.is_point() ) {
      in_conflict = is_interior_in_conflict_none_sp(p, q, r, s, t, tag);
    } else { // both p and q are segments
      in_conflict = is_interior_in_conflict_none_ss(p, q, r, s, t, tag);
    }

    return in_conflict;
  }

  //------------------------------------------------------------------------


  bool
  is_interior_in_conflict_none_pp(const Site_2& p, const Site_2& q,
				  const Site_2& , const Site_2& ,
				  const Site_2& t, Method_tag ) const
  {
    CGAL_precondition( p.is_point() && q.is_point() && t.is_point() );
    CGAL_USE(p);
    CGAL_USE(q);
    CGAL_USE(t);
    return false;
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none_ps(const Site_2& sp, const Site_2& sq,
				  const Site_2& r, const Site_2& s,
				  const Site_2& st, Method_tag ) const
  {
    CGAL_precondition( sp.is_point() && sq.is_segment() && st.is_point() );

    if ( same_points(sp, sq.source_site()) ||
	 same_points(sp, sq.target_site()) ) {
      return false;
    }

    Line_2 lq = compute_supporting_line(sq.supporting_site());

    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    Point_2 p = sp.point(), t = st.point();

    Line_2 lperp = compute_linf_perpendicular(lq, t);

    Oriented_side op = oriented_side_of_line(lq, p);
    Oriented_side ot = oriented_side_of_line(lq, t);

    bool on_same_side =
      ((op == ON_POSITIVE_SIDE && ot == ON_POSITIVE_SIDE) ||
       (op == ON_NEGATIVE_SIDE && ot == ON_NEGATIVE_SIDE));

    Comparison_result res =
      compare_linf_distances_to_line(lq, t, p);

    Oriented_side opqr = vpqr.oriented_side(lperp);
    Oriented_side oqps = vqps.oriented_side(lperp);

    bool on_different_side =
      ((opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE) ||
       (opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE));

    return ( on_same_side && (res == SMALLER) && on_different_side );
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none_sp(const Site_2& p, const Site_2& q,
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag tag) const
  {
    return is_interior_in_conflict_none_ps(q, p, s, r, t, tag);
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none_ss(const Site_2& p, const Site_2& q,
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag ) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() && t.is_point() );

    Voronoi_vertex_2 vpqr(p, q, r);
    Voronoi_vertex_2 vqps(q, p, s);

    Line_2 lp = compute_supporting_line(p.supporting_site());
    Line_2 lq = compute_supporting_line(q.supporting_site());

    // first orient lp according to its Voronoi vertices
    if ( vpqr.is_degenerate_Voronoi_circle() ) {
      Site_2 tpqr = Site_2::construct_site_2(vpqr.degenerate_point());

      if ( same_points(tpqr, p.source_site()) ||
	   same_points(tpqr, p.target_site()) ) {
	if ( vqps.oriented_side(lp) != ON_POSITIVE_SIDE ) {
	  lp = opposite_line(lp);
	}
      }
    } else {
      if ( vpqr.oriented_side(lp) != ON_POSITIVE_SIDE ) {
	lp = opposite_line(lp);
      }
    }
#if 0 // OLD CODE
    if (  ( vpqr.is_degenerate_Voronoi_circle() &&
	    same_points(vpqr.degenerate_point(), p.source_site()) ) ||
	  ( vpqr.is_degenerate_Voronoi_circle() &&
	    same_points(vpqr.degenerate_point(), p.target_site()) )  ) {
      if ( vqps.oriented_side(lp) != ON_POSITIVE_SIDE ) {
	lp = opposite_line(lp);
      }
    } else {
      if ( vpqr.oriented_side(lp) != ON_POSITIVE_SIDE ) {
	lp = opposite_line(lp);
      }
    }
#endif

    // then orient lq according to its Voronoi vertices
    if ( vpqr.is_degenerate_Voronoi_circle() ) {
      Site_2 tpqr = Site_2::construct_site_2(vpqr.degenerate_point());

      if ( same_points(tpqr, q.source_site()) ||
	   same_points(tpqr, q.target_site()) ) {
	if ( vqps.oriented_side(lq) != ON_POSITIVE_SIDE ) {
	  lq = opposite_line(lq);
	}
      }
    } else {
      if ( vpqr.oriented_side(lq) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    }
#if 0 // OLD CODE
    if (  ( vpqr.is_degenerate_Voronoi_circle() &&
	    same_points(vpqr.degenerate_point(), q.source_site()) ) ||
	  ( vpqr.is_degenerate_Voronoi_circle() &&
	    same_points(vpqr.degenerate_point(), q.target_site()) )  ) {
      if ( vqps.oriented_side(lq) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    } else {
      if ( vpqr.oriented_side(lq) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    }
#endif

    Point_2 tp = t.point();

    // check if t is on the same side as the Voronoi vertices
    Oriented_side ot_lp = oriented_side_of_line(lp, tp);
    Oriented_side ot_lq = oriented_side_of_line(lq, tp);

    if ( ot_lp != ON_POSITIVE_SIDE || ot_lq != ON_POSITIVE_SIDE ) {
      return false;
    }

    Line_2 lperp;

    Comparison_result res =
      compare_linf_distances_to_lines(tp, lp, lq);

    if ( res == SMALLER ) {
      lperp = compute_linf_perpendicular(lp, tp);
    } else {
      lperp = compute_linf_perpendicular(lq, tp);
    }

    CGAL_precondition( ot_lp != ON_ORIENTED_BOUNDARY &&
		       ot_lq != ON_ORIENTED_BOUNDARY );

    // check of lperp separates the two Voronoi vertices
    Oriented_side opqr_perp = vpqr.oriented_side(lperp);
    Oriented_side oqps_perp = vqps.oriented_side(lperp);

    bool on_different_side =
      (opqr_perp == ON_POSITIVE_SIDE &&
       oqps_perp == ON_NEGATIVE_SIDE) ||
      (opqr_perp == ON_NEGATIVE_SIDE &&
       oqps_perp == ON_POSITIVE_SIDE);

    return ( on_different_side );
  }

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

public:
  typedef Boolean           result_type;
  typedef Site_2            argument_type;

  Boolean   operator()(const Site_2& p, const Site_2& q, const Site_2& r,
		       const Site_2& s, const Site_2& t, Sign sgn) const
  {
    Boolean retval;
    CGAL_SDG_DEBUG(
        std::cout << "debug finite-edge-int-cf entering (p,q,r,s,t,sgn)= "
        << p << ' ' << q << ' ' << r << ' ' << s << ' ' << t
        << ' ' << sgn << std::endl; );
    if ( sgn == POSITIVE ) {
      retval = is_interior_in_conflict_none(p, q, r, s, t, Method_tag());
    } else if ( sgn == NEGATIVE ) {
      retval = is_interior_in_conflict_both(p, q, r, s, t, Method_tag());
    } else {
      retval = is_interior_in_conflict_touch(p, q, r, s, t, Method_tag());
    }
    CGAL_SDG_DEBUG( std::cout <<
        "debug finite-edge-int-cf with (p,q,r,s,t,sgn)= "
        << p << ' ' << q << ' ' << r << ' ' << s << ' ' << t
        << ' ' << sgn << " returns " << retval << std::endl; );
    return retval;
  }


  // philaris: tocheck
  Boolean   operator()(const Site_2& p, const Site_2& q, const Site_2& r,
		       const Site_2& t, Sign sgn) const
  {
    CGAL_SDG_DEBUG(
        std::cout << "debug finite-edge-int-cf entering (p,q,r,t,sgn)= ("
        << p << ") (" << q << ") (" << r <<  " (r not ignored)) ("
        << t << ")  "
        << sgn << std::endl; );

    /*
    if ( t.is_point() ) {
      return ( sgn == NEGATIVE );
    }
    */

    // if p, q, t are points, then return (sgn == NEGATIVE)
    if (t.is_point() && p.is_point() && q.is_point()) {
      CGAL_SDG_DEBUG( std::cout <<
          "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
          << p << ") (" << q << ") (" << r <<  " ) ("
          << t << ")  "
          << sgn << " returns " << (sgn==NEGATIVE) << std::endl; );
      return ( sgn == NEGATIVE );
    }

    if (t.is_point() && (sgn == NEGATIVE)) {
      bool is_p_pnt = p.is_point();
      bool is_q_pnt = q.is_point();
      bool parabola_case = (is_p_pnt && (! is_q_pnt)) ||
                           (is_q_pnt && (! is_p_pnt))   ;
      if (parabola_case) {
        if (is_p_pnt) {
          Point_2 pp = p.point();
          Line_2 lq = compute_supporting_line(q.supporting_site());
          if ( oriented_side_of_line(lq, pp)
              == ON_NEGATIVE_SIDE ) {
            lq = opposite_line(lq);
          }
          Voronoi_vertex_2 vpqr(p, q, r);
          Line_2 lqperp = compute_linf_perpendicular(lq, p.point());
          Oriented_side opqr = vpqr.oriented_side(lqperp);
          Boolean on_different_parabola_arcs = (opqr == ON_NEGATIVE_SIDE);

          //if ( !on_different_parabola_arcs ) { return true; }
          if (certainly( ! on_different_parabola_arcs ) ) {
            CGAL_SDG_DEBUG( std::cout
                << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << true << std::endl; );
            return true;
          }
          if (! is_certain( ! on_different_parabola_arcs ) ) {
            return indeterminate<Boolean>();
          }
        }
        if (is_q_pnt) {
          Point_2 qq = q.point();
          Line_2 lp = compute_supporting_line(p.supporting_site());
          if ( oriented_side_of_line(lp, qq)
              == ON_NEGATIVE_SIDE ) {
            lp = opposite_line(lp);
          }
          Line_2 lpperp = compute_linf_perpendicular(lp, qq);
          Voronoi_vertex_2 vpqr(p, q, r);
          Oriented_side opqr = vpqr.oriented_side(lpperp);
          Boolean on_different_parabola_arcs =
            (opqr == ON_POSITIVE_SIDE);

          //if ( !on_different_parabola_arcs ) { return true; }
          if (certainly( ! on_different_parabola_arcs ) ) {
            CGAL_SDG_DEBUG( std::cout
                << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << true << std::endl; );
            return true;
          }
          if (! is_certain( ! on_different_parabola_arcs ) ) {
            return indeterminate<Boolean>();
          }
        }
        CGAL_SDG_DEBUG( std::cout <<
            "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
            << p << ") (" << q << ") (" << r <<  ") ("
            << t << ")  "
            << sgn << " returns " << false << std::endl; );
        return false;
      }

      CGAL_SDG_DEBUG( std::cout <<
          "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
          << p << ") (" << q << ") (" << r <<  ") ("
          << t << ")  "
          << sgn << " returns " << true << std::endl; );
      return true;
    }

    // here and on, t is a segment
    // or t is a point and sgn != NEGATIVE

    if ( sgn != NEGATIVE ) {
      // philaris: tocheck

      // philaris: here, I assert that sgn is never zero
      // tocheck
      CGAL_assertion( sgn == POSITIVE );

      CGAL_assertion( p.is_point() || q.is_point() );
      CGAL_assertion( p.is_segment() || q.is_segment() );

      CGAL_assertion(t.is_point());

      if ( p.is_point() && q.is_segment() ) {
        if (same_points(p, q.source_site()) ||
            same_points(p, q.target_site())   ) {
          CGAL_SDG_DEBUG( std::cout <<
              "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
              << p << ") (" << q << ") (" << r <<  ") ("
              << t << ")  "
              << sgn << " returns " << false << std::endl; );
          return false;
        } else { // p is not endpoint of q

          CGAL_SDG_DEBUG( std::cout << "debug fecf p is not endpoint of q"
            << std::endl; );

          CGAL_assertion( ! is_site_h_or_v(q) ) ;

          Line_2 lseg = compute_supporting_line(q.supporting_site());

          Oriented_side ossegp =
            oriented_side_of_line(lseg, p.point());

          Oriented_side ossegt =
            oriented_side_of_line(lseg, t.point());

          CGAL_assertion(ossegp != ON_ORIENTED_BOUNDARY);
          CGAL_assertion(ossegt != ON_ORIENTED_BOUNDARY);

          if (ossegp != ossegt) {
            CGAL_SDG_DEBUG( std::cout <<
                "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << false << std::endl; );
            return false;
          } else {
            // here, t and p are on the same side of q

            CGAL_SDG_DEBUG(
                std::cout << "debug fecf: t and p on the same side of q"
                << std::endl; );

            // compute linf projection of infinite
            // Voronoi vertex (p, inf, q)

            Point_2 proj_of_infv;
            bool has_lseg_neg_slope =
              CGAL::sign(lseg.a()) == CGAL::sign(lseg.b());

            if (has_lseg_neg_slope) {
              proj_of_infv =
                compute_horizontal_projection(lseg, p.point());
            } else {
              proj_of_infv =
                compute_vertical_projection(lseg, p.point());
            }

            // compute linf projection of pqr Voronoi vertex
            Voronoi_vertex_2 vpqr(p, q, r);

            Point_2 proj_of_vpqr =
              compute_linf_projection_nonhom(lseg, vpqr.point());

            Line_2 line_inf =
              compute_linf_perpendicular(lseg, proj_of_infv);

            Line_2 line_fin =
              compute_linf_perpendicular(lseg, proj_of_vpqr);

            Oriented_side os_line_inf =
              oriented_side_of_line(line_inf, t.point());

            Oriented_side os_line_fin =
              oriented_side_of_line(line_fin, t.point());

            Boolean retval = ( os_line_inf == os_line_fin) ? false : true;

            CGAL_SDG_DEBUG( std::cout <<
                "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << retval << std::endl; );

            return retval;

          } // end of case where
            // t, p are on the same side of q
        } // end of case where p is not endpoint of q
      } // end of case where p is point and q is segment

      if ( p.is_segment() && q.is_point() ) {
        if (same_points(q, p.source_site()) ||
            same_points(q, p.target_site())   ) {
          CGAL_SDG_DEBUG( std::cout <<
              "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
              << p << ") (" << q << ") (" << r <<  ") ("
              << t << ")  "
              << sgn << " returns " << false << std::endl; );
          return false;
        } else { // q is not endpoint of p
          CGAL_assertion( ! is_site_h_or_v(p) );

          Line_2 lseg = compute_supporting_line(p.supporting_site());

          Oriented_side ossegq =
            oriented_side_of_line(lseg, q.point());

          Oriented_side ossegt =
            oriented_side_of_line(lseg, t.point());

          CGAL_assertion(ossegq != ON_ORIENTED_BOUNDARY);
          CGAL_assertion(ossegt != ON_ORIENTED_BOUNDARY);

          if (ossegq != ossegt) {
            CGAL_SDG_DEBUG( std::cout <<
              "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
              << p << ") (" << q << ") (" << r <<  ") ("
              << t << ")  "
              << sgn << " returns " << false << std::endl; );
            return false;
          } else {
            // here, t and q are on the same side of p

            // compute linf projection of infinite
            // Voronoi vertex (p, inf, q)

            Point_2 proj_of_infv;
            bool has_lseg_neg_slope =
              CGAL::sign(lseg.a()) == CGAL::sign(lseg.b());

            if (has_lseg_neg_slope) {
              proj_of_infv =
                compute_vertical_projection(lseg, q.point());
            } else {
              proj_of_infv =
                compute_horizontal_projection(lseg, q.point());
            }

            // compute linf projection of pqr Voronoi vertex
            Voronoi_vertex_2 vpqr(p, q, r);

            Point_2 proj_of_vpqr =
              compute_linf_projection_nonhom(lseg, vpqr.point());

            Line_2 line_inf =
              compute_linf_perpendicular(lseg, proj_of_infv);

            Line_2 line_fin =
              compute_linf_perpendicular(lseg, proj_of_vpqr);

            Oriented_side os_line_inf =
              oriented_side_of_line(line_inf, t.point());

            Oriented_side os_line_fin =
              oriented_side_of_line(line_fin, t.point());

            Boolean retval = ( os_line_inf == os_line_fin) ? false : true;

            CGAL_SDG_DEBUG( std::cout <<
                "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << retval << std::endl; );

            return retval;

          } // end of case where
            // t, q are on the same side of p
        } // end of case where q is not endpoint of p
      } // end of case where q is point and p is segment

      CGAL_SDG_DEBUG( std::cout <<
          "debug finite-edge-int-cf with (p,q,r,t,sgn)= ("
          << p << ") (" << q << ") (" << r <<  ") ("
          << t << ")  "
          << sgn << " returns " << false << std::endl; );

      return false;
    }

    // here and on, sgn == NEGATIVE

    if ( p.is_segment() || q.is_segment() ) {
      // philaris: p and q are consecutive on convex hull
      // therefore, one of them is a point and the other is a segment
      CGAL_assertion( p.is_point() || q.is_point() );

      // philaris:

      CGAL_assertion( p.is_point() ?
                      ((! same_points(p, q.source_site())) &&
                       (! same_points(p, q.target_site())))    :
                      ((! same_points(q, p.source_site())) &&
                       (! same_points(q, p.target_site())))
                    ) ;

      //sandeep: if t is segment and it not intersects the wedge
      //as defined by s and p in the basic predicate
      //intersects_segment_interior_inf_wedge_sp
      // then return false
      if (t.is_segment()) {
        bool result;

        if (p.is_point()) {
          CGAL_assertion(q.is_segment());
          //p may be end point of t
          if ( same_points(p,t.source_site())
               || same_points(p,t.target_site()) ) {
            CGAL_SDG_DEBUG(
                std::cout
                << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << true << std::endl; );
            return true;
          } else {
            if (same_points(q.source_site(), t.source_site()) ||
                same_points(q.source_site(), t.target_site()) ||
                same_points(q.target_site(), t.source_site()) ||
                same_points(q.target_site(), t.target_site())   ) {
              result = true;
            } else {
              Point_2 pp = p.point();
              Line_2 lq = compute_supporting_line(q.supporting_site());
              if ( oriented_side_of_line(lq, pp)
                   == ON_NEGATIVE_SIDE ) {
                lq = opposite_line(lq);
              }
              Line_2 lqperp = compute_linf_perpendicular(lq, pp);
              Voronoi_vertex_2 vpqr(p, q, r);
              Oriented_side opqr = vpqr.oriented_side(lqperp);
              Boolean on_different_parabola_arcs =
                (opqr == ON_NEGATIVE_SIDE);

              //if ( !on_different_parabola_arcs ) { return true; }
              if (certainly( ! on_different_parabola_arcs ) ) {
                CGAL_SDG_DEBUG(
                    std::cout
                    << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
                    << p << ") (" << q << ") (" << r <<  ") ("
                    << t << ")  "
                    << sgn << " returns " << true << std::endl; );
                return true;
              }
              if (! is_certain( ! on_different_parabola_arcs ) ) {
                return indeterminate<Boolean>();
              }

              // here, we are on different parabola arcs

              CGAL_SDG_DEBUG( std::cout
                  << "debug fecf pqrt on diff parab arcs" << std::endl; );

              CGAL_SDG_DEBUG( std::cout
                 << "debug fecf call intersects( " << q << ", "
                 << p << ", " << t << ")" << std::endl; ) ;
              result = intersects_segment_interior_inf_wedge_sp(q,p,t);
            }
          }
        } else { // p is segment and q is point
          CGAL_assertion(q.is_point());
          //q may be endpoint of t
          if (   same_points(q,t.source_site())
              || same_points(q,t.target_site())) {
            CGAL_SDG_DEBUG( std::cout
                << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
                << p << ") (" << q << ") (" << r <<  ") ("
                << t << ")  "
                << sgn << " returns " << true << std::endl; );
            return true;
          } else {
            if (same_points(p.source_site(), t.source_site()) ||
                same_points(p.source_site(), t.target_site()) ||
                same_points(p.target_site(), t.source_site()) ||
                same_points(p.target_site(), t.target_site())   ) {
              result = true;
            } else {
              Point_2 qq = q.point();
              Line_2 lp = compute_supporting_line(p.supporting_site());
              if ( oriented_side_of_line(lp, qq)
                   == ON_NEGATIVE_SIDE ) {
                lp = opposite_line(lp);
              }
              Line_2 lpperp = compute_linf_perpendicular(lp, qq);
              Voronoi_vertex_2 vpqr(p, q, r);
              Oriented_side opqr = vpqr.oriented_side(lpperp);
              Boolean on_different_parabola_arcs =
                (opqr == ON_POSITIVE_SIDE);

              //if ( !on_different_parabola_arcs ) { return true; }
              if (certainly( ! on_different_parabola_arcs ) ) {
                CGAL_SDG_DEBUG(
                    std::cout
                    << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
                    << p << ") (" << q << ") (" << r <<  ") ("
                    << t << ")  "
                    << sgn << " returns " << true << std::endl; );
                return true;
              }
              if (! is_certain( ! on_different_parabola_arcs ) ) {
                return indeterminate<Boolean>();
              }

              // here, we are on different parabola arcs

              CGAL_SDG_DEBUG( std::cout
                  << "debug fecf pqrt on diff parab arcs" << std::endl; );

              CGAL_SDG_DEBUG( std::cout
                 << "debug fecf call intersects( " << p << ", "
                 << q << ", " << t << ")" << std::endl; ) ;
              result = intersects_segment_interior_inf_wedge_sp(p,q,t);
            }
          }
        }

        if (result == false) {
          CGAL_SDG_DEBUG( std::cout
              << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
              << p << ") (" << q << ") (" << r <<  ") ("
              << t << ")  "
              << sgn << " returns " << false << std::endl; );
          return false;
        }
      }//end of t is segment

      // philaris: tocheck
      CGAL_SDG_DEBUG( std::cout
          << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
          << p << ") (" << q << ") (" << r <<  ") ("
          << t << ")  "
          << sgn << " returns " << true << std::endl; );
      // philaris: always return true
      return true;
    }

    // here p and q are points and s is a segment

    bool p_is_endpoint =
      same_points(p, t.source_site()) || same_points(p, t.target_site());
    bool q_is_endpoint =
      same_points(q, t.source_site()) || same_points(q, t.target_site());

    // philaris: change for Linf
    if ((cmpx(p.point(), q.point()) == EQUAL) ||
        (cmpy(p.point(), q.point()) == EQUAL)    ) {
      // if p or q have one same coordinate, behave like in L2
      return ( p_is_endpoint && q_is_endpoint );
    } else {
      // philaris: tocheck
      CGAL_assertion(sgn == NEGATIVE);
      CGAL_SDG_DEBUG( std::cout
          << "debug finite-edge-int-cf tocheck (p,q,r,t,sgn)= ("
          << p << ") (" << q << ") (" << r <<  ") ("
          << t << ")  "
          << sgn << " returns " << true << std::endl; );
      return true;
    }
  }

  Boolean   operator()(const Site_2& p, const Site_2& q, const Site_2& t,
		       Sign sgn) const
  {
    CGAL_SDG_DEBUG(
        std::cout << "debug finite-edge-int-cf entering (p,q,t,sgn)= "
        << p << ' ' << q << ' ' << t << ' '
        << "(sgn " << sgn << " not ignored)"  << std::endl; );

    CGAL_assertion( ! ( p.is_segment() && q.is_segment()) );

    if ( p.is_segment() || q.is_segment() ) {
      Site_2 siteseg = (p.is_point())? q: p;

      if ( is_site_h_or_v(siteseg) ) {
        if (p.is_point()) {
          CGAL_assertion( same_points(p, q.source_site()) ||
                          same_points(p, q.target_site())   ) ;
        } else {
          CGAL_assertion( same_points(q, p.source_site()) ||
                          same_points(q, p.target_site())   ) ;
        }
        CGAL_SDG_DEBUG(
            std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
            << p << ' ' << q << ' ' << t << ' '
            << "(sgn " << sgn << " not ignored)" << " returns "
            << false << std::endl; );
        return false;
      }

      CGAL_assertion(! is_site_h_or_v(siteseg));
      // t is segment
      if (t.is_segment()) {
        CGAL_assertion(sgn == NEGATIVE);
        // sandeep if t does not intersect the wedge as defined by
        //intersects_segment_interior_inf_wedge_sp(seg,pt,t) then return true
        // this works if p or q is point
        bool result;
        if (p.is_point()) { // p is point, q is segment
          CGAL_assertion((! same_points(p, q.source_site())) &&
                         (! same_points(p, q.target_site()))    );
          result = intersects_segment_interior_inf_wedge_sp(q,p,t);
        } else { // q is point and p is segment
          CGAL_assertion((! same_points(q, p.source_site())) &&
                         (! same_points(q, p.target_site()))    );
          result = intersects_segment_interior_inf_wedge_sp(p,q,t);
        }

        CGAL_SDG_DEBUG(
            std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
            << p << ' ' << q << ' ' << t << ' '
            << "(sgn " << sgn << " not ignored)" << " returns "
            << result << std::endl; );

        if (result == true) {
          return true;
        } else {
          return false;
        }
        //return false;
      }

      CGAL_assertion(t.is_point());

      if (sgn == POSITIVE) {
        if (p.is_point() && q.is_segment()) {
          if (same_points(p, q.source_site()) ||
              same_points(p, q.target_site())   ) {
            CGAL_SDG_DEBUG(
                std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
                << p << ' ' << q << ' ' << t << ' '
                << "(sgn " << sgn << " not ignored)" << " returns "
                << false << std::endl; );
            return false;
          }
        } // end of case: p: point, q: segment
        else
        // if (p.is_segment() && q.is_point())
        {
          if (same_points(q, p.source_site()) ||
              same_points(q, p.target_site())   ) {
            CGAL_SDG_DEBUG(
                std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
                << p << ' ' << q << ' ' << t << ' '
                << "(sgn " << sgn << " not ignored)" << " returns "
                << false << std::endl; );
            return false;
          }
        } // end of case: p: segment, q: point

        Point_2 pnt = (p.is_point())? p.point(): q.point();
        Line_2 lseg = (p.is_point()) ?
          compute_supporting_line(q.supporting_site()) :
          compute_supporting_line(p.supporting_site()) ;

        Oriented_side ossegpnt =
          oriented_side_of_line(lseg, pnt);

        Oriented_side ossegt =
          oriented_side_of_line(lseg, t.point());

        CGAL_assertion(ossegpnt != ON_ORIENTED_BOUNDARY);
        CGAL_assertion(ossegt != ON_ORIENTED_BOUNDARY);

        if (ossegpnt != ossegt) {
          CGAL_SDG_DEBUG(
              std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
              << p << ' ' << q << ' ' << t << ' '
              << "(sgn " << sgn << " not ignored)" << " returns "
              << false << std::endl; );
          return false;
        } else {
          // here, t and pnt are on the same side of seg

          Point_2 proj_hor =
            compute_horizontal_projection(lseg, pnt);

          Point_2 proj_ver =
            compute_vertical_projection(lseg, pnt);

          Line_2 lineh =
            compute_linf_perpendicular(lseg, proj_hor);

          Line_2 linev =
            compute_linf_perpendicular(lseg, proj_ver);

          Oriented_side os_lineh =
            oriented_side_of_line(lineh, t.point());

          Oriented_side os_linev =
            oriented_side_of_line(linev, t.point());

          Boolean retval = ( os_lineh == os_linev ) ? false : true;

          CGAL_SDG_DEBUG(
              std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
              << p << ' ' << q << ' ' << t << ' '
              << "(sgn " << sgn << " not ignored)" << " returns "
              << retval << std::endl; );

          return retval;

        } // end of case: t and pnt on the same side of seg

      } // end of case: sgn == POSITIVE

      return false;
    } // end of case: any of p, q is a segment

    CGAL_SDG_DEBUG(
        std::cout << "debug tocheck here p, q points" << std::endl; );

    // both p and q are points
    if ( t.is_point() ) {
      // philaris: adaptation as follows:
      // predicate: position of t related to bbox(p, q)

      Comparison_result cmpxpt = cmpx(p.point(), t.point());
      Comparison_result cmpxtq = cmpx(t.point(), q.point());
      Comparison_result cmpypt = cmpy(p.point(), t.point());
      Comparison_result cmpytq = cmpy(t.point(), q.point());

      //Sign s1 = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
      Sign s1 = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

      CGAL_assertion( s1 != ZERO );

      CGAL_SDG_DEBUG( std::cout << "debug finite-edge-int-cf (p,q,t)= "
        << p << ' ' << q << ' ' << t
        << "  s1= " << s1 << std::endl; );

      CGAL_SDG_DEBUG(
          std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
          << p << ' ' << q << ' ' << t << ' '
          << "(sgn " << sgn << " not ignored)" << " returns "
          << ( s1 == NEGATIVE ) << std::endl; );

      return ( s1 == NEGATIVE );
    }

    CGAL_SDG_DEBUG( std::cout << "debug tocheck" << std::endl; );

    bool bp =
      same_points(p, t.source_site()) || same_points(p, t.target_site());
    bool bq =
      same_points(q, t.source_site()) || same_points(q, t.target_site());

    Comparison_result cmpxpq = cmpx(p.point(), q.point());
    Comparison_result cmpypq = cmpy(p.point(), q.point());

    if ((cmpxpq == EQUAL) || (cmpypq == EQUAL)) {
      CGAL_SDG_DEBUG(
          std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
          << p << ' ' << q << ' ' << t << ' '
          << "(sgn " << sgn << " not ignored)" << " returns "
          << ( bp && bq ) << std::endl; );
      return ( bp && bq );
    } else {
      // check if interior of segment has non-empty intersection
      // with the bounding box of p, q
      Boolean retval = intersects_segment_interior_bbox(t, p, q);
      CGAL_SDG_DEBUG(
          std::cout << "debug finite-edge-int-cf with (p,q,t,sgn)= "
          << p << ' ' << q << ' ' << t << ' '
          << "(sgn " << sgn << " not ignored)" << " returns "
          << retval << std::endl; );
      return retval;
    }
  }

};


//------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_FINITE_EDGE_INTERIOR_CONFLICT_C2_H

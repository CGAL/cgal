// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>




#ifndef CGAL_SVD_FINITE_EDGE_INTERIOR_2_H
#define CGAL_SVD_FINITE_EDGE_INTERIOR_2_H

#include <CGAL/predicates/Svd_basic_predicates_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>
#include <CGAL/predicates/Svd_are_same_points_C2.h>
#include <CGAL/predicates/Svd_are_same_segments_C2.h>

CGAL_BEGIN_NAMESPACE


//-----------------------------------------------------------------------------

template<class K, class Method_tag>
class Svd_finite_edge_interior_2
  : public Svd_basic_predicates_C2<K>
{
public:

  typedef Svd_basic_predicates_C2<K>          Base;
  typedef Svd_voronoi_vertex_2<K,Method_tag>  Voronoi_vertex_2;
  
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Line_2               Line_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;

  typedef typename Base::Homogeneous_point_2  Homogeneous_point_2;

private:
  typedef Svd_are_same_points_C2<K>           Are_same_points_2;
  typedef Svd_are_same_segments_C2<K>         Are_same_segments_2;

  typedef typename K::Intersections_tag       ITag;

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

private:

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Site_2& p, const Site_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag tag) const
  {
    bool in_conflict(false);

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
				  const Site_2& t, Method_tag tag) const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() );

    Point_2 p = sp.point(), q = sq.point();

    if ( t.is_point() ) { return true; }

    Line_2 lt = compute_supporting_line(t.supporting_segment());

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
      compare_squared_distances_to_line(lt, p, q);

    if ( res == EQUAL ) { return true; }

    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    
    Line_2 lperp;
    if ( res == SMALLER ) {
      // p is closer to lt than q
      lperp = compute_perpendicular(lt, p);
    } else {
      // q is closer to lt than p
      lperp = compute_perpendicular(lt, q);
    }

    Oriented_side opqr = vpqr.oriented_side(lperp);
    Oriented_side oqps = vqps.oriented_side(lperp);

    return ( opqr == oqps );
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both_ss(const Site_2& p, const Site_2& q,
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() );
    return true;
  }

  //--------------------------------------------------------------------

  bool
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

  bool
  is_interior_in_conflict_both_ps_p(const Site_2& p, const Site_2& q,
				    const Site_2& r, const Site_2& s,
				    const Site_2& t, Method_tag tag) const
  {
    CGAL_precondition( t.is_point() );

    //    Line_2 lq = compute_supporting_line(q);
    Line_2 lq = compute_supporting_line(q.supporting_segment());

    Comparison_result res =
      compare_squared_distances_to_line(lq, p.point(), t.point());

    if ( res != SMALLER ) { return true; }

    Voronoi_vertex_2 vpqr(p, q, r);
    Voronoi_vertex_2 vqps(q, p, s);

    Line_2 lperp = compute_perpendicular(lq, p.point());
      
    Oriented_side opqr = vpqr.oriented_side(lperp);
    Oriented_side oqps = vqps.oriented_side(lperp);

    return (opqr == oqps);
  }

  //--------------------------------------------------------------------

  bool check_if_exact(const Site_2& t1, const Tag_false&) const
  {
    return true;
  }

  bool check_if_exact(const Site_2& t1, const Tag_true&) const
  {
    return t1.is_exact();
  }

  bool
  is_interior_in_conflict_both_ps_s(const Site_2& sp, const Site_2& sq,
				    const Site_2& r, const Site_2& s,
				    const Site_2& st, Method_tag tag) const
  {
    CGAL_precondition( st.is_segment() );
    Point_2 p = sp.point();
    //    Segment_2 q = sq.segment(), t = st.segment();

    Line_2 lq = compute_supporting_line(sq.supporting_segment());

    if ( oriented_side_of_line(lq, p) == ON_NEGATIVE_SIDE ) {
      lq = opposite_line(lq); 
    }

    if ( same_points(sp, st.source_site()) ||
	 same_points(sp, st.target_site()) ) {

      Line_2 lqperp = compute_perpendicular(lq, p);

      Voronoi_vertex_2 vpqr(sp, sq, r);
      Voronoi_vertex_2 vqps(sq, sp, s);

      Oriented_side opqr = vpqr.oriented_side(lqperp);
      Oriented_side oqps = vqps.oriented_side(lqperp);

      bool on_different_parabola_arcs =
	((opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE) ||
	 (opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE));

      if ( !on_different_parabola_arcs ) { return true; }

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
	compare_squared_distances_to_line(lq, p, t1.point());

      return ( res == LARGER );
    }

    Line_2 lt = compute_supporting_line(st.supporting_segment());

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

    Line_2 l = compute_perpendicular(lt, p);

    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    Oriented_side o_l_pqr = vpqr.oriented_side(l);
    Oriented_side o_l_qps = vqps.oriented_side(l);
    if ( o_l_pqr == ON_POSITIVE_SIDE &&
	 o_l_qps == ON_NEGATIVE_SIDE ) { return false; }
    if ( o_l_pqr == ON_NEGATIVE_SIDE &&
	 o_l_qps == ON_POSITIVE_SIDE ) { return true; }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>> HERE I NEED TO CHECK THE BOUNDARY CASES <<<<<<
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line_2 lqperp = compute_perpendicular(lq, p);

    Oriented_side opqr = vpqr.oriented_side(lqperp);
    Oriented_side oqps = vqps.oriented_side(lqperp);

    bool on_different_parabola_arcs =
      ((opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE) ||
       (opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE));

    if ( !on_different_parabola_arcs ) { return true; }
      
    Homogeneous_point_2 pv = projection_on_line(lq, p);
    Homogeneous_point_2 hp(p);

    pv = midpoint(pv, hp);

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

  bool
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

    if ( (p.is_point() && q.is_point()) ||
	 (p.is_segment() && q.is_segment()) ) { 
      return true;
    }

    if ( p.is_point() && q.is_segment() ) {
      Line_2 lq = compute_supporting_line(q.supporting_segment());
    
      Comparison_result res =
	compare_squared_distances_to_line(lq, p.point(), t.point());

      return (res != SMALLER);
    }

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
				  const Site_2& r, const Site_2& s,
				  const Site_2& t, Method_tag tag) const
  {
    CGAL_precondition( p.is_point() && q.is_point() && t.is_point() );
    return false;
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none_ps(const Site_2& sp, const Site_2& sq,
				  const Site_2& r, const Site_2& s,
				  const Site_2& st, Method_tag tag) const
  {
    CGAL_precondition( sp.is_point() && sq.is_segment() && st.is_point() );

    if ( same_points(sp, sq.source_site()) ||
	 same_points(sp, sq.target_site()) ) {
      return false;
    }
   
    Line_2 lq = compute_supporting_line(sq.supporting_segment());

    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    Point_2 p = sp.point(), t = st.point();

    Line_2 lperp = compute_perpendicular(lq, t);

    Oriented_side op = oriented_side_of_line(lq, p);
    Oriented_side ot = oriented_side_of_line(lq, t);

    bool on_same_side =
      ((op == ON_POSITIVE_SIDE && ot == ON_POSITIVE_SIDE) ||
       (op == ON_NEGATIVE_SIDE && ot == ON_NEGATIVE_SIDE));

    Comparison_result res =
      compare_squared_distances_to_line(lq, t, p);

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
				  const Site_2& t, Method_tag tag) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() && t.is_point() );

    Voronoi_vertex_2 vpqr(p, q, r);
    Voronoi_vertex_2 vqps(q, p, s);

    Line_2 lp = compute_supporting_line(p.supporting_segment());
    Line_2 lq = compute_supporting_line(q.supporting_segment());

    // first orient lp according to its Voronoi vertices
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

    // then orient lq according to its Voronoi vertices
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

    Point_2 tp = t.point();

    // check if t is on the same side as the Voronoi vertices
    Oriented_side ot_lp = oriented_side_of_line(lp, tp);
    Oriented_side ot_lq = oriented_side_of_line(lq, tp);

    if ( ot_lp != ON_POSITIVE_SIDE || ot_lq != ON_POSITIVE_SIDE ) {
      return false;
    }

    Line_2 lperp;

    Comparison_result res =
      compare_squared_distances_to_lines(tp, lp, lq);

    if ( res == SMALLER ) {
      lperp = compute_perpendicular(lp, tp);
    } else {
      lperp = compute_perpendicular(lq, tp);
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
  bool operator()(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& s, const Site_2& t, Sign sgn) const
  {
    if ( sgn == POSITIVE ) {
      return is_interior_in_conflict_none(p, q, r, s, t, Method_tag());
    } else if ( sgn == NEGATIVE ) {
      return is_interior_in_conflict_both(p, q, r, s, t, Method_tag());
    } else {
      return is_interior_in_conflict_touch(p, q, r, s, t, Method_tag());
    }
  }


  bool operator()(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, Sign sgn) const
  {
    if ( t.is_point() ) {
      return ( sgn == NEGATIVE );
    }

    if ( sgn != NEGATIVE ) {
      return false;
    }

    if ( p.is_segment() || q.is_segment() ) {
      return false;
    }

    bool p_is_endpoint =
      same_points(p, t.source_site()) || same_points(p, t.target_site());
    bool q_is_endpoint =
      same_points(q, t.source_site()) || same_points(q, t.target_site());

    return ( p_is_endpoint && q_is_endpoint );
  }

  bool operator()(const Site_2& p, const Site_2& q, const Site_2& t,
		  Sign sgn) const
  {
    if ( p.is_segment() || q.is_segment()) {
      return false;
    }

    // both p and q are points
    if ( t.is_point() ) {
      RT dtpx = p.point().x() - t.point().x();
      RT minus_dtpy = -p.point().y() + t.point().y();
      RT dtqx = q.point().x() - t.point().x();
      RT dtqy = q.point().y() - t.point().y();

      Sign s1 = sign_of_determinant2x2(dtpx, minus_dtpy, dtqy, dtqx);

      CGAL_assertion( s1 != ZERO );
      return ( s1 == NEGATIVE );
    }

    bool bp =
      same_points(p, t.source_site()) || same_points(p, t.target_site());
    bool bq =
      same_points(q, t.source_site()) || same_points(q, t.target_site());
						       

    return ( bp && bq );
  }

};


//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_FINITE_EDGE_INTERIOR_2_H

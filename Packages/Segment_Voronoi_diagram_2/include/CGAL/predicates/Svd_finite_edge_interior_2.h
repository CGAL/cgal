#ifndef CGAL_SVD_FINITE_EDGE_INTERIOR_2_H
#define CGAL_SVD_FINITE_EDGE_INTERIOR_2_H

#include <CGAL/predicates/Svd_basic_predicates_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>


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
      in_conflict =      
	is_interior_in_conflict_both(p.point(), q.point(),
				     r, s, t, tag);
    } else if ( p.is_segment() && q.is_segment() ) {
      in_conflict =
	is_interior_in_conflict_both(p.segment(), q.segment(),
				     r, s, t, tag);
    } else if ( p.is_point() && q.is_segment() ) {
      in_conflict =
	is_interior_in_conflict_both(p.point(), q.segment(),
				     r, s, t, tag);
    } else { // p is a segment and q is a point
      in_conflict =
	is_interior_in_conflict_both(p.segment(), q.point(),
				     r, s, t, tag);
    }

    return in_conflict;
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Point_2& p, const Point_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag tag) const
  {
    if ( t.is_point() ) { return true; }

    Line_2 lt = compute_supporting_line(t.segment());

    Oriented_side op, oq;

    if ( p == t.segment().source() || 
	 p == t.segment().target() ) {
      op = ON_ORIENTED_BOUNDARY;
    } else {
      op = oriented_side_of_line(lt, p);
    }

    if ( q == t.segment().source() || 
	 q == t.segment().target() ) {
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

    Site_2 sp(p), sq(q);
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
  is_interior_in_conflict_both(const Segment_2& p, const Segment_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag) const
  {
    return true;
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Point_2& p, const Segment_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag tag) const
  {
    if ( p == q.source() || p == q.target() ) {
      return false;
    }   

    if ( t.is_point() ) {
      return is_interior_in_conflict_both(p, q, r, s, t.point(), tag);
    }
    return is_interior_in_conflict_both(p, q, r, s, t.segment(), tag);
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Point_2& p, const Segment_2& q,
			       const Site_2& r, const Site_2& s,
			       const Point_2& t, Method_tag tag) const
  {
    Line_2 lq = compute_supporting_line(q);

    Comparison_result res =
      compare_squared_distances_to_line(lq, p, t);

    if ( res != SMALLER ) { return true; }

    Site_2 sp(p), sq(q);
    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    Line_2 lperp = compute_perpendicular(lq, p);
      
    Oriented_side opqr = vpqr.oriented_side(lperp);
    Oriented_side oqps = vqps.oriented_side(lperp);

    return (opqr == oqps);
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Point_2& p, const Segment_2& q,
			       const Site_2& r, const Site_2& s,
			       const Segment_2& t, Method_tag tag) const
  {
    Line_2 lt = compute_supporting_line(t);
    Line_2 lq = compute_supporting_line(q);

    if ( oriented_side_of_line(lq, p) == ON_NEGATIVE_SIDE ) {
      lq = opposite_line(lq); 
    }

    if ( p == t.source() || p == t.target() ) {
      Line_2 lqperp = compute_perpendicular(lq, p);


      Site_2 sp(p), sq(q);
      Voronoi_vertex_2 vpqr(sp, sq, r);
      Voronoi_vertex_2 vqps(sq, sp, s);

      Oriented_side opqr = vpqr.oriented_side(lqperp);
      Oriented_side oqps = vqps.oriented_side(lqperp);

      bool on_different_parabola_arcs =
	((opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE) ||
	 (opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE));

      if ( !on_different_parabola_arcs ) { return true; }

      Point_2 t1;
      if ( p == t.source() ) {
	t1 = t.target();
      } else {
	t1 = t.source();
      }

      Oriented_side o_p = oriented_side_of_line(lq, p);
      Oriented_side o_t1 = oriented_side_of_line(lq, t1);

      if ( (o_p == ON_POSITIVE_SIDE && o_t1 == ON_NEGATIVE_SIDE) ||
	   (o_p == ON_NEGATIVE_SIDE && o_t1 == ON_POSITIVE_SIDE) ) {
	return true;
      }
	     

      Comparison_result res =
	compare_squared_distances_to_line(lq, p, t1);

      return ( res == LARGER );
    }

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

    Site_2 sp(p), sq(q);
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
  is_interior_in_conflict_both(const Segment_2& p, const Point_2& q,
			       const Site_2& r, const Site_2& s,
			       const Site_2& t, Method_tag tag) const
  {
    return is_interior_in_conflict_both(q, p, s, r, t, tag);
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
      Line_2 lq = compute_supporting_line(q.segment());
    
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

    Point_2 tp = t.point();

    bool in_conflict(false);

    if ( p.is_point() && q.is_point() ) {
      in_conflict =
	is_interior_in_conflict_none(p.point(), q.point(),
				     r, s, tp, tag);
    } else if ( p.is_point() && q.is_segment() ) {
      in_conflict =
	is_interior_in_conflict_none(p.point(), q.segment(),
				     r, s, tp, tag);
    } else if ( p.is_segment() && q.is_point() ) {
      in_conflict =
	is_interior_in_conflict_none(p.segment(), q.point(),
				     r, s, tp, tag);
    } else { // both p and q are segments
      in_conflict =
	is_interior_in_conflict_none(p.segment(), q.segment(),
				     r, s, tp, tag);
    }

    return in_conflict;
  }

  //------------------------------------------------------------------------


  bool
  is_interior_in_conflict_none(const Point_2& p, const Point_2& q,
			       const Site_2& r, const Site_2& s,
			       const Point_2& t, Method_tag tag) const
  {
    return false;
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none(const Point_2& p, const Segment_2& q,
			       const Site_2& r, const Site_2& s,
			       const Point_2& t, Method_tag tag) const
  {
    if ( p == q.source() || p == q.target() ) {
      return false;
    }
   
    Line_2 lq = compute_supporting_line(q);

    Site_2 sp(p), sq(q);
    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

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
  is_interior_in_conflict_none(const Segment_2& p, const Point_2& q,
			       const Site_2& r, const Site_2& s,
			       const Point_2& t, Method_tag tag) const
  {
    return is_interior_in_conflict_none(q, p, s, r, t, tag);
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none(const Segment_2& p, const Segment_2& q,
			       const Site_2& r, const Site_2& s,
			       const Point_2& t, Method_tag tag) const
  {
    Site_2 sp(p), sq(q);
    Voronoi_vertex_2 vpqr(sp, sq, r);
    Voronoi_vertex_2 vqps(sq, sp, s);

    Line_2 lp = compute_supporting_line(p);
    Line_2 lq = compute_supporting_line(q);

    // first orient lp according to its Voronoi vertices
    if (  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.degenerate_point() == p.source() ) ||
	  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.degenerate_point() == p.target() )  ) {
      //      CGAL_assertion
      //	( !vqps.is_same_point(p.source(), tag) &&
      //	  !vqps.is_same_point(p.target(), tag) );
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
	    vpqr.degenerate_point() == q.source() ) ||
	  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.degenerate_point() == q.target() )  ) {
      //      CGAL_assertion
      //	( !vqps.is_same_point(q.source(), tag) &&
      //	  !vqps.is_same_point(q.target(), tag) );
      if ( vqps.oriented_side(lq) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    } else {
      if ( vpqr.oriented_side(lq) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    }

    // check if t is on the same side as the Voronoi vertices
    Oriented_side ot_lp = oriented_side_of_line(lp, t);
    Oriented_side ot_lq = oriented_side_of_line(lq, t);

    if ( ot_lp != ON_POSITIVE_SIDE || ot_lq != ON_POSITIVE_SIDE ) {
      return false;
    }

    Line_2 lperp;

    Comparison_result res =
      compare_squared_distances_to_lines(t, lp, lq);

    if ( res == SMALLER ) {
      lperp = compute_perpendicular(lp, t);
    } else {
      lperp = compute_perpendicular(lq, t);
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
    bool res;
    if ( sgn == POSITIVE ) {
      res = is_interior_in_conflict_none(p, q, r, s, t, Method_tag());
    } else if ( sgn == NEGATIVE ) {
      res = is_interior_in_conflict_both(p, q, r, s, t, Method_tag());
    } else {
      res = is_interior_in_conflict_touch(p, q, r, s, t, Method_tag());
    }

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
      (p.point() == t.source() || p.point() == t.target());
    bool q_is_endpoint =
      (q.point() == t.source() || q.point() == t.target());

    return ( p_is_endpoint && q_is_endpoint );
  }

  bool operator()(const Site_2& p, const Site_2& q, const Site_2& t,
		  Sign sgn) const
  {
    if ( sgn != ZERO ) {
      return false;
    }

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

    bool bp = ( (p.point() == t.source()) || (p.point() == t.target()) );
    bool bq = ( (q.point() == t.source()) || (q.point() == t.target()) );

    return ( bp && bq );
  }

};


//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_FINITE_EDGE_INTERIOR_2_H

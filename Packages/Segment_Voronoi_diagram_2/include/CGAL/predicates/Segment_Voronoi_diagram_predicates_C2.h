#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H


#include <CGAL/Segment_Voronoi_diagram_constructions_C2.h>
#include <CGAL/Segment_Voronoi_diagram_site_2.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< class K >
class Svd_are_same_points_2
{
public:
  typedef typename K::Point_2      Point_2;

  typedef typename K::Compare_x_2  compare_x_2;
  typedef typename K::Compare_y_2  compare_y_2;
public:

  inline
  bool operator()(const Point_2& p, const Point_2& q) const
  {
    if ( compare_x_2()(p, q) != EQUAL ) { return false; }
    Comparison_result res = compare_y_2()(p, q);
    return res == EQUAL;
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< class K >
class Svd_is_endpoint_of_segment_2
{
public:
  typedef typename K::Point_2                       Point_2;
  typedef typename K::Segment_2                     Segment_2;
  typedef typename CGAL::Svd_are_same_points_2<K>   Are_same_points_2;

  typedef typename K::Compare_x_2  compare_x_2;
  typedef typename K::Compare_y_2  compare_y_2;
public:

  inline
  bool operator()(const Point_2& p, const Segment_2& s) const
  {
    Are_same_points_2 are_same_points;
    return ( are_same_points(p, s.source()) ||
	     are_same_points(p, s.target()) );
  }
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


template<class R, class Method_tag>
class Svd_incircle_2
{
private:
  typedef typename R::Point_2            Point;
  typedef typename R::Segment_2          Segment;
  typedef typename R::Site_2             Site;
  typedef CGAL::Svd_voronoi_vertex_2<R>  Voronoi_vertex;

  typedef typename R::FT                 FT;
  typedef typename R::RT                 RT;

private:
  Sign incircle(const Site& p, const Site& q,
		const Point& t) const
  {
    if ( p.is_point() && q.is_point() ) {
      Orientation o = orientation(p.point(), q.point(), t);

      if ( o != COLLINEAR ) {
	return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
      }

      
      RT dtpx = p.point().x() - t.x();
      RT dtpy = p.point().y() - t.y();
      RT dtqx = q.point().x() - t.x();
      RT minus_dtqy = -q.point().y() + t.y();
      
      Sign s = sign_of_determinant2x2(dtpx, dtpy, minus_dtqy, dtqx);

      CGAL_assertion( s != ZERO );

      return s;
    }

    CGAL_assertion( p.is_point() || q.is_point() );

    Orientation o;
    if ( p.is_point() && q.is_segment() ) {
      Point pq = (p.point() == q.source()) ? q.target() : q.source();
      o = orientation(p.point(), pq, t);
    } else { // p is a segment and q is a point
      Point pp = (q.point() == p.source()) ? p.target() : p.source();
      o = orientation(pp, q.point(), t);
    }
    return ( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE;
  }

  //-----------------------------------------------------------------------


  Sign incircle(const Point& p, const Point& q,
		const Segment& t) const
  {    
    if ( (p == t.source() || p == t.target()) &&
	 (q == t.source() || q == t.target()) ) {
	// if t is the segment joining p and q then t must be a vertex
	// on the convex hull
	return NEGATIVE;
    } else if ( p == t.source() || p == t.target() ) {
      // p is an endpoint of t
      // in this case the p,q,oo vertex is destroyed only if the
      // other endpoint of t is beyond
      Point pt = (p == t.source()) ? t.target() : t.source();
      Orientation o = orientation(p, q, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else if ( q == t.source() || q == t.target() ) {
      Point pt = (q == t.source()) ? t.target() : t.source();
      Orientation o = orientation(p, q, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // maybe here I should immediately return POSITIVE;
      // since we insert endpoints of segments first, p and q cannot
      // be consecutive points on the convex hull if one of the
      // endpoints of t is to the right of the line pq.
      Orientation o1 = orientation(p, q, t.source());
      Orientation o2 = orientation(p, q, t.target());

      if ( o1 == RIGHT_TURN || o2 == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return POSITIVE;
    }
  }


  Sign incircle(const Segment& p, const Point& q,
		const Segment& t) const
  {
    if ( q == t.source() && q == t.target() ) {
      Point pp = (q == p.source()) ? p.target() : p.source();
      Point pt = (q == t.source()) ? t.target() : t.source();

      Orientation o = orientation(pp, q, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      return POSITIVE;
    }
  }


  Sign incircle(const Point& p, const Segment& q,
		const Segment& t) const
  {
    if ( p == t.source() || p == t.target() ) {
      Point pq = (p == q.source()) ? q.target() : q.source();
      Point pt = (p == t.source()) ? t.target() : t.source();

      Orientation o = orientation(p, pq, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // if p is not an endpoint of t, then either p and q should
      // not be on the convex hull or t does not affect the vertex
      // of p and q.
      return POSITIVE;
    }
  }


  Sign incircle(const Site& p, const Site& q,
		const Segment& t) const
  {
    if ( p.is_point() && q.is_point() ) {
      return incircle(p.point(), q.point(), t);
    } else if ( p.is_point() && q.is_segment() ) {
      return incircle(p.point(), q.segment(), t);
    } else { // p is a segment and q is a point
      return incircle(p.segment(), q.point(), t);
    }
  }


public:

  Sign operator()(const Site& p, const Site& q,
		  const Site& r, const Site& t) const
  {
    Voronoi_vertex v(p, q, r);

    return v.incircle(t, Method_tag());
  }


  

  inline Sign operator()(const Site& p, const Site& q,
			 const Site& t) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    if ( p.is_point() && q.is_segment() ) {
      // p must be an endpoint of q
      CGAL_assertion( (p.point() == q.source()) ||
		      (p.point() == q.target()) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( (p.source() == q.point()) ||
		      (p.target() == q.point()) );
    }

    if ( t.is_point() ) {
      return incircle(p, q, t.point());
    }

    return incircle(p, q, t.segment());
  }


};



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------



template<class K, class Method_tag>
class Svd_finite_edge_interior_2
{
public:
  typedef typename K::Point_2            Point;
  typedef typename K::Segment_2          Segment;
  typedef typename K::Line_2             Line;
  typedef CGAL::Svd_voronoi_vertex_2<K>  Voronoi_vertex;
  typedef typename K::Site_2             Site;
  typedef typename K::FT                 FT;
  typedef typename K::RT                 RT;

  typedef typename Voronoi_vertex::vertex_t   vertex_t;
  typedef typename Voronoi_vertex::PPP_Type   PPP_Type;
  typedef typename Voronoi_vertex::PPS_Type   PPS_Type;
  typedef typename Voronoi_vertex::PSS_Type   PSS_Type;
  typedef typename Voronoi_vertex::SSS_Type   SSS_Type;

  typedef typename Voronoi_vertex::Sqrt_1     Sqrt_1;
  typedef typename Voronoi_vertex::Sqrt_2     Sqrt_2;
  typedef typename Voronoi_vertex::Sqrt_3     Sqrt_3;

#if 0
  class Line
  {
  private:
    RT a_, b_, c_;

  public:
    Line() : a_(1), b_(0), c_(0) {}
    Line(const RT& a, const RT& b, const RT& c)
      : a_(a), b_(b), c_(c) {}

    RT a() const { return a_; }
    RT b() const { return b_; }
    RT c() const { return c_; }

  };
#endif


  struct Homogeneous_point_2
  {
    RT hx_, hy_, hw_;

    Homogeneous_point_2() : hx_(0), hy_(0), hw_(1) {}
    Homogeneous_point_2(const RT& hx, const RT& hy, const RT& hw)
      :  hx_(hx), hy_(hy), hw_(hw)
    {
      CGAL_precondition( !(CGAL_NTS is_zero(hw_)) );
    }

    Homogeneous_point_2(const Point& p)
      : hx_(p.x()), hy_(p.y()), hw_(1) {}

    Homogeneous_point_2(const Homogeneous_point_2& other)
      : hx_(other.hx_), hy_(other.hy_), hw_(other.hw_) {}

    RT hx() const { return hx_; }
    RT hy() const { return hy_; }
    RT hw() const { return hw_; }

    FT x() const { return hx_ / hw_; }
    FT y() const { return hy_ / hw_; }
  };

public:
  static
  Homogeneous_point_2
  projection_on_line(const Line& l, const Point& p)
  {
    RT ab = l.a() * l.b();

    RT hx = CGAL_NTS square(l.b()) * p.x()
      - ab * p.y() - l.a() * l.c();
    RT hy = CGAL_NTS square(l.a()) * p.y()
      - ab * p.x() - l.b() * l.c();
    RT hw = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());

    return Homogeneous_point_2(hx, hy, hw);
  }


  static
  Homogeneous_point_2
  midpoint(const Point& p1, const Point& p2)
  {
    RT hx = p1.x() + p2.x();
    RT hy = p1.y() + p2.y();
    RT hw = RT(2);

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Homogeneous_point_2
  midpoint(const Homogeneous_point_2& p1,
	   const Homogeneous_point_2& p2)
  {
    RT hx = p1.hx() * p2.hw() + p2.hx() * p1.hw();
    RT hy = p1.hy() * p2.hw() + p2.hy() * p1.hw();
    RT hw = RT(2) * p1.hw() * p2.hw();

    return Homogeneous_point_2(hx, hy, hw);
  }
  

  static
  Line compute_supporting_line(const Segment& s)
  {
#if 1
    RT a, b, c;
    Voronoi_vertex::compute_supporting_line(s, a, b, c);
    return Line(a, b, c);
#else
    return Line(s);
#endif
  }

  static
  Comparison_result
  compare_squared_distances_to_line(const Line& l, const Point& p,
				    const Point& q)
  {
#if 1
    RT d2_lp = CGAL_NTS square(l.a() * p.x() + l.b() * p.y() + l.c());
    RT d2_lq = CGAL_NTS square(l.a() * q.x() + l.b() * q.y() + l.c());

    return CGAL_NTS compare(d2_lp, d2_lq);
#else
    FT d2_lp = squared_distance_to_line(p, l);
    FT d2_lq = squared_distance_to_line(q, l);

    return CGAL_NTS compare(d2_lp, d2_lq);
#endif
  }

  static
  Comparison_result
  compare_squared_distances_to_lines(const Point& p, const Line& l1,
				     const Line& l2)
  {
#if 1
    RT d2_l1 = CGAL_NTS square(l1.a() * p.x() + l1.b() * p.y() + l1.c());
    RT d2_l2 = CGAL_NTS square(l2.a() * p.x() + l2.b() * p.y() + l2.c());

    RT n1 = CGAL_NTS square(l1.a()) + CGAL_NTS square(l1.b());
    RT n2 = CGAL_NTS square(l2.a()) + CGAL_NTS square(l2.b());

    return CGAL_NTS compare(d2_l1 * n2, d2_l2 * n1);
#else

    FT d2_p_from_l1 = squared_distance_to_line(p, l1);
    FT d2_p_from_l2 = squared_distance_to_line(p, l2);

    return CGAL_NTS compare(d2_p_from_l1, d2_p_from_l2);
#endif
  }

  static
  Oriented_side
  oriented_side_of_line(const Line& l, const Point& p)
  {
#if 1
    Sign s = CGAL_NTS sign(l.a() * p.x() + l.b() * p.y() + l.c());
    if ( s == ZERO ) { return ON_ORIENTED_BOUNDARY; }
    return ( s == POSITIVE ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
#else
    return l.oriented_side(p);
#endif
  }


  static
  Oriented_side
  oriented_side_of_line(const Line& l, const Homogeneous_point_2& p)
  {
    Sign s1 =
      CGAL_NTS sign(l.a() * p.hx() + l.b() * p.hy() + l.c() * p.hw());
    Sign s_hw = CGAL_NTS sign(p.hw());

    Sign s = Sign(s1 * s_hw);

    if ( s == ZERO ) { return ON_ORIENTED_BOUNDARY; }
    return ( s == POSITIVE ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  }



  static
  Line compute_perpendicular(const Line& l, const Point& p)
  {
#if 1
    RT a, b, c;
    a = -l.b();
    b = l.a();
    c = l.b() * p.x() - l.a() * p.y();
    return Line(a, b, c);
#else
    return l.perpendicular(p);
#endif
  }

  static
  Line opposite_line(const Line& l)
  {
#if 1
    return Line(-l.a(), -l.b(), -l.c());
#else
    return l.opposite();
#endif
  }

private:
  static
  FT squared_distance_to_line(const Point& p, const Line& l)
  {
    Point proj = l.projection(p);
    Segment s(p, proj);
    return s.squared_length();
  }

  static
  FT
  squared_distance_to_line(const Point& p, const Segment& s)
  {
    return squared_distance_to_line(p, Line(s));
  }


private:

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
#if 1
  bool
  is_interior_in_conflict_both(const Site& p, const Site& q,
			       const Site& r, const Site& s,
			       const Site& t, Method_tag tag) const
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

#if 0
    CGAL::Object o = make_object(tag);

    Svd_finite_edge_interior_2<K,Sqrt_field_tag>  P;

    Sqrt_field_tag tag1;
    if ( CGAL::assign(tag1, o) ) {
    } else {
      bool in_conflict2 = P(p, q, r, s, t, NEGATIVE);
      CGAL_assertion( in_conflict == in_conflict2 );
    }
#endif

    return in_conflict;
  }
#endif
  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Point& p, const Point& q,
			       const Site& r, const Site& s,
			       const Site& t, Method_tag tag) const
  {
    if ( t.is_point() ) { return true; }

    Line lt = compute_supporting_line(t.segment());

    // MK:: THIS IS REALLY IMPORTANT BECAUSE IT DEPENDS ON HOW LINES
    //    ARE IMPLEMENTED *******************************************
    //      change from passing a line to passing a triple of numbers
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    Site sp(p), sq(q);
    Voronoi_vertex vpqr(sp, sq, r);
    Voronoi_vertex vqps(sq, sp, s);

    
    Line lperp;
    if ( res == SMALLER ) {
      // p is closer to lt than q
      lperp = compute_perpendicular(lt, p);
    } else {
      // q is closer to lt than p
      lperp = compute_perpendicular(lt, q);
    }

    Oriented_side opqr = vpqr.oriented_side(lperp, tag);
    Oriented_side oqps = vqps.oriented_side(lperp, tag);

    return ( opqr == oqps );
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Segment& p, const Segment& q,
			       const Site& r, const Site& s,
			       const Site& t, Method_tag) const
  {
    return true;
  }

  //--------------------------------------------------------------------

  bool
  is_interior_in_conflict_both(const Point& p, const Segment& q,
			       const Site& r, const Site& s,
			       const Site& t, Method_tag tag) const
  {
    if ( p == q.source() || p == q.target() ) {
      return false;
    }   

    if ( t.is_point() ) {
      return is_interior_in_conflict_both(p, q, r, s, t.point(), tag);
    }
    return is_interior_in_conflict_both(p, q, r, s, t.segment(), tag);
  }

  bool
  is_interior_in_conflict_both(const Point& p, const Segment& q,
			       const Site& r, const Site& s,
			       const Point& t, Method_tag tag) const
  {
    Line lq = compute_supporting_line(q);

    Comparison_result res =
      compare_squared_distances_to_line(lq, p, t);

    if ( res != SMALLER ) { return true; }

    Site sp(p), sq(q);
    Voronoi_vertex vpqr(sp, sq, r);
    Voronoi_vertex vqps(sq, sp, s);

    Line lperp = compute_perpendicular(lq, p);
      
    Oriented_side opqr = vpqr.oriented_side(lperp, tag);
    Oriented_side oqps = vqps.oriented_side(lperp, tag);

    return (opqr == oqps);
  }

  bool
  is_interior_in_conflict_both(const Point& p, const Segment& q,
			       const Site& r, const Site& s,
			       const Segment& t, Method_tag tag) const
  {
    Line lt = compute_supporting_line(t);
    Line lq = compute_supporting_line(q);

    if ( oriented_side_of_line(lq, p) == ON_NEGATIVE_SIDE ) {
      lq = opposite_line(lq); 
    }

    if ( p == t.source() || p == t.target() ) {
      Line lqperp = compute_perpendicular(lq, p);


      Site sp(p), sq(q);
      Voronoi_vertex vpqr(sp, sq, r);
      Voronoi_vertex vqps(sq, sp, s);

      Oriented_side opqr = vpqr.oriented_side(lqperp, tag);
      Oriented_side oqps = vqps.oriented_side(lqperp, tag);

      bool on_different_parabola_arcs =
	((opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE) ||
	 (opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE));

      if ( !on_different_parabola_arcs ) { return true; }

      Point t1;
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
      CGAL_NTS compare(lt.a() * lq.b(), lt.b() * lq.a());
    bool are_parallel = (res == EQUAL);
      
    if ( are_parallel ) {
      Sign sgn = CGAL_NTS sign(lt.a() * lq.a() + lt.b() * lq.b());
      bool have_opposite_directions = (sgn == NEGATIVE);
      if ( have_opposite_directions ) { lq = opposite_line(lq); }

      if ( oriented_side_of_line(lq, p) == oriented_side_of_line(lt, p) ) {
	return true;
      }

      if ( have_opposite_directions ) {
	lq = opposite_line(lq); 	  
      }
    }

    Line l = compute_perpendicular(lt, p);

    Site sp(p), sq(q);
    Voronoi_vertex vpqr(sp, sq, r);
    Voronoi_vertex vqps(sq, sp, s);

    Oriented_side o_l_pqr = vpqr.oriented_side(l, tag);
    Oriented_side o_l_qps = vqps.oriented_side(l, tag);
    if ( o_l_pqr == ON_POSITIVE_SIDE &&
	 o_l_qps == ON_NEGATIVE_SIDE ) { return false; }
    if ( o_l_pqr == ON_NEGATIVE_SIDE &&
	 o_l_qps == ON_POSITIVE_SIDE ) { return true; }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>> HERE I NEED TO CHECK THE BOUNDARY CASES <<<<<<
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    Line lqperp = compute_perpendicular(lq, p);

    Oriented_side opqr = vpqr.oriented_side(lqperp, tag);
    Oriented_side oqps = vqps.oriented_side(lqperp, tag);

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
  is_interior_in_conflict_both(const Segment& p, const Point& q,
			       const Site& r, const Site& s,
			       const Site& t, Method_tag tag) const
  {
    return is_interior_in_conflict_both(q, p, s, r, t, tag);
  }


  //--------------------------------------------------------------------

#if 0
  bool is_interior_in_conflict_both(const Site& p, const Site& q,
				    const Site& r, const Site& s,
				    const Site& t, Ring_tag tag)
    const
  {
    // THIS NEEDS TO BE CORRECTED
    return
      is_interior_in_conflict_both(p, q, r, s, t, Sqrt_field_tag());
  }

  bool is_interior_in_conflict_both(const Site& p, const Site& q,
				    const Site& r, const Site& s,
				    const Site& t, Sqrt_field_tag tag)
    const
  {
    // checks if interior of voronoi edge is in conflict if both extrema 
    // of the voronoi edge are in conclift
    // return true if interior is in conflict; false otherwise

    if ( p.is_point() && q.is_point() ) {
      // both p and q are points
      if ( t.is_point() ) { return true; }

      Line lt = compute_supporting_line(t.segment());

      Oriented_side op, oq;

      if ( p.point() == t.segment().source() || 
	   p.point() == t.segment().target() ) {
	op = ON_ORIENTED_BOUNDARY;
      } else {
	op = oriented_side_of_line(lt, p.point());
      }

      if ( q.point() == t.segment().source() || 
	   q.point() == t.segment().target() ) {
	oq = ON_ORIENTED_BOUNDARY;
      } else {
	oq = oriented_side_of_line(lt, q.point());
      }
    

      if ((op == ON_POSITIVE_SIDE && oq == ON_NEGATIVE_SIDE) ||
	  (op == ON_NEGATIVE_SIDE && oq == ON_POSITIVE_SIDE) ||
	  (op == ON_ORIENTED_BOUNDARY || oq == ON_ORIENTED_BOUNDARY)) {
	return true;
      }

      Comparison_result res =
	compare_squared_distances_to_line(lt, p.point(), q.point());
      //      FT d2_p_from_lt = squared_distance_to_line(p.point(), lt);
      //      FT d2_q_from_lt = squared_distance_to_line(q.point(), lt);

      //      res = CGAL_NTS compare(d2_p_from_lt, d2_q_from_lt);

      if ( res == EQUAL ) { return true; }

      Voronoi_vertex vpqr(p, q, r);
      Voronoi_vertex vqps(q, p, s);

      Line lperp;
      if ( res == SMALLER ) {
	// p is closer to lt than q
	lperp = compute_perpendicular(lt, p.point());
      } else {
	// q is closer to lt than p
	lperp = compute_perpendicular(lt, q.point());
      }

      Oriented_side opqr = vpqr.oriented_side(lperp, tag);
      Oriented_side oqps = vqps.oriented_side(lperp, tag);

      return ( opqr == oqps );
    } else if ( p.is_segment() && q.is_segment() ) {
      // both p and q are segments
      return true;

    } else if ( p.is_point() && q.is_segment() ) {
      // p is a point q is a segment
      if ( p.point() == q.segment().source() || 
	   p.point() == q.segment().target() ) {
	return false;
      }    
      if ( t.is_point() ) {
	Line lq = compute_supporting_line(q.segment());

	//	FT d2_p_from_lq = squared_distance_to_line(p.point(), lq);
	//	FT d2_t_from_lq = squared_distance_to_line(t.point(), lq);

	//	Comparison_result res =
	//	  CGAL_NTS compare(d2_p_from_lq, d2_t_from_lq);

	Comparison_result res =
	  compare_squared_distances_to_line(lq, p.point(), t.point());

	if ( res != SMALLER ) { return true; }

	Voronoi_vertex vpqr(p, q, r);
	Voronoi_vertex vqps(q, p, s);

	//	Line lperp = lq.perpendicular(p.point());
	Line lperp = compute_perpendicular(lq, p.point());

	//	Oriented_side opqr = lperp.oriented_side(vpqr.point());
	//	Oriented_side oqps = lperp.oriented_side(vqps.point());

	Oriented_side opqr = vpqr.oriented_side(lperp, tag);
	Oriented_side oqps = vqps.oriented_side(lperp, tag);

	return (opqr == oqps);
      } else {
	// t is a segment
	Line lt = compute_supporting_line(t.segment());
	Line lq = compute_supporting_line(q.segment());

	//	if ( lq.has_on_negative_side(p.point()) ) { 
	if ( oriented_side_of_line(lq, p.point()) ==
	     ON_NEGATIVE_SIDE ) {
	  lq = opposite_line(lq);
	  //	  lq = lq.opposite(); 
	}

	if ( p.point() == t.source() || p.point() == t.target() ) {
	  //	  Line lqperp = lq.perpendicular(p.point());
	  Line lqperp = compute_perpendicular(lq, p.point());

	  Voronoi_vertex vpqr(p, q, r);
	  Voronoi_vertex vqps(q, p, s);

	  //	  Oriented_side opqr = lqperp.oriented_side(vpqr.point());
	  //	  Oriented_side oqps = lqperp.oriented_side(vqps.point());

	  Oriented_side opqr = vpqr.oriented_side(lqperp, tag);
	  Oriented_side oqps = vqps.oriented_side(lqperp, tag);

	  bool on_different_parabola_arcs =
	    ((opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE) ||
	     (opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE));

	  if ( !on_different_parabola_arcs ) { return true; }

	  Point t1;
	  if ( p.point() == t.source() ) {
	    t1 = t.target();
	  } else {
	    t1 = t.source();
	  }

	  //	  Oriented_side o_p = lq.oriented_side(p.point());
	  //	  Oriented_side o_t1 = lq.oriented_side(t1);

	  Oriented_side o_p = oriented_side_of_line(lq, p.point());
	  Oriented_side o_t1 = oriented_side_of_line(lq, t1);

	  if ( (o_p == ON_POSITIVE_SIDE && o_t1 == ON_NEGATIVE_SIDE) ||
	       (o_p == ON_NEGATIVE_SIDE && o_t1 == ON_POSITIVE_SIDE) ) {
	    return true;
	  }
	     
	  //	  FT d2_p_from_lq = squared_distance_to_line(p.point(), lq);
	  //	  FT d2_t1_from_lq = squared_distance_to_line(t1, lq);

	  //	  Comparison_result res =
	  //	    CGAL_NTS compare(d2_p_from_lq, d2_t1_from_lq);

	  Comparison_result res =
	    compare_squared_distances_to_line(lq, p.point(), t1);

	  return ( res == LARGER );
	}

	//	if ( lt.has_on_negative_side(p.point()) ) {
	if ( oriented_side_of_line(lt, p.point()) == 
	     ON_NEGATIVE_SIDE) {
	  //	  lt = lt.opposite();
	  lt = opposite_line(lt);
	}


	Comparison_result res =
	  CGAL_NTS compare(lt.a() * lq.b(), lt.b() * lq.a());
	bool are_parallel = (res == EQUAL);
      
	if ( are_parallel ) {
	  Sign sgn = CGAL_NTS sign(lt.a() * lq.a() + lt.b() * lq.b());
	  bool have_opposite_directions = (sgn == NEGATIVE);
	  if ( have_opposite_directions ) {
	    //  lq = lq.opposite();
	    lq = opposite_line(lq);
	  }

	  //	  if ( lq.oriented_side(p.point()) == 
	  //	       lt.oriented_side(p.point()) ) {
	  if ( oriented_side_of_line(lq, p.point()) == 
	       oriented_side_of_line(lt, p.point()) ) {
	    return true;
	  }
	  if ( have_opposite_directions ) {
	    lq = opposite_line(lq);
	    //	    lq = lq.opposite(); 	  
	  }
	}

	//	Line l = lt.perpendicular(p.point());
	Line l = compute_perpendicular(lt, p.point());

	Voronoi_vertex vpqr(p, q, r);
	Voronoi_vertex vqps(q, p, s);

	//	Oriented_side o_l_pqr = l.oriented_side(vpqr.point());
	//	Oriented_side o_l_qps = l.oriented_side(vqps.point());

	Oriented_side o_l_pqr = vpqr.oriented_side(l, tag);
	Oriented_side o_l_qps = vqps.oriented_side(l, tag);


	if ( o_l_pqr == ON_POSITIVE_SIDE &&
	     o_l_qps == ON_NEGATIVE_SIDE ) { return false; }
	if ( o_l_pqr == ON_NEGATIVE_SIDE &&
	     o_l_qps == ON_POSITIVE_SIDE ) { return true; }

	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//>>>>>>>>>> HERE I NEED TO CHECK THE BOUNDARY CASES <<<<<<
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	//	Line lqperp = lq.perpendicular(p.point());
	Line lqperp = compute_perpendicular(lq, p.point());

	//	Oriented_side opqr = lqperp.oriented_side(vpqr.point());
	//	Oriented_side oqps = lqperp.oriented_side(vqps.point());

	Oriented_side opqr = vpqr.oriented_side(lqperp, tag);
	Oriented_side oqps = vqps.oriented_side(lqperp, tag);

	bool on_different_parabola_arcs =
	  ((opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE) ||
	   (opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE));

	if ( !on_different_parabola_arcs ) { return true; }
      

	Point pv = lq.projection(p.point());
	pv = midpoint(pv, p.point());

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
    } else {
      // q is a point p is a segment
      return is_interior_in_conflict_both(q, p, s, r, t, tag);
    }
  }

#endif

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_touch(const Site& p, const Site& q,
				const Site& r, const Site& s,
				const Site& t, Method_tag tag) const
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
      Line lq = compute_supporting_line(q.segment());
    
      Comparison_result res =
	compare_squared_distances_to_line(lq, p.point(), t.point());

      return (res != SMALLER);
    }

    return is_interior_in_conflict_touch(q, p, s, r, t, tag);
  }

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------


#if 1

  bool
  is_interior_in_conflict_none(const Site& p, const Site& q,
			       const Site& r, const Site& s,
			       const Site& t, Method_tag tag) const
  {
    if ( t.is_segment() ) { return false; }

    Point tp = t.point();

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
#endif

  //------------------------------------------------------------------------


  bool
  is_interior_in_conflict_none(const Point& p, const Point& q,
			       const Site& r, const Site& s,
			       const Point& t, Method_tag tag) const
  {
    return false;
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none(const Point& p, const Segment& q,
			       const Site& r, const Site& s,
			       const Point& t, Method_tag tag) const
  {
    if ( p == q.source() || p == q.target() ) {
      return false;
    }
   
    Line lq = compute_supporting_line(q);

    Site sp(p), sq(q);
    Voronoi_vertex vpqr(sp, sq, r);
    Voronoi_vertex vqps(sq, sp, s);

    Line lperp = compute_perpendicular(lq, t);

    Oriented_side op = oriented_side_of_line(lq, p);
    Oriented_side ot = oriented_side_of_line(lq, t);

    bool on_same_side =
      ((op == ON_POSITIVE_SIDE && ot == ON_POSITIVE_SIDE) ||
       (op == ON_NEGATIVE_SIDE && ot == ON_NEGATIVE_SIDE));

    Comparison_result res =
      compare_squared_distances_to_line(lq, t, p);

    Oriented_side opqr = vpqr.oriented_side(lperp, tag);
    Oriented_side oqps = vqps.oriented_side(lperp, tag);

    bool on_different_side =
      ((opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE) ||
       (opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE));

    return ( on_same_side && (res == SMALLER) && on_different_side );
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none(const Segment& p, const Point& q,
			       const Site& r, const Site& s,
			       const Point& t, Method_tag tag) const
  {
    return is_interior_in_conflict_none(q, p, s, r, t, tag);
  }

  //------------------------------------------------------------------------

  bool
  is_interior_in_conflict_none(const Segment& p, const Segment& q,
			       const Site& r, const Site& s,
			       const Point& t, Method_tag tag) const
  {
    Site sp(p), sq(q);
    Voronoi_vertex vpqr(sp, sq, r);
    Voronoi_vertex vqps(sq, sp, s);

    Line lp = compute_supporting_line(p);
    Line lq = compute_supporting_line(q);

    // first orient lp according to its Voronoi vertices
    if (  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.point() == p.source() ) ||
	  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.point() == p.target() )  ) {
      //      CGAL_assertion
      //	( !vqps.is_same_point(p.source(), tag) &&
      //	  !vqps.is_same_point(p.target(), tag) );
      if ( vqps.oriented_side(lp, tag) != ON_POSITIVE_SIDE ) {
	lp = opposite_line(lp);
      }
    } else {
      if ( vpqr.oriented_side(lp, tag) != ON_POSITIVE_SIDE ) {
	lp = opposite_line(lp);
      }
    }

    // then orient lq according to its Voronoi vertices
    if (  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.point() == q.source() ) ||
	  ( vpqr.is_degenerate_Voronoi_circle() &&
	    vpqr.point() == q.target() )  ) {
      //      CGAL_assertion
      //	( !vqps.is_same_point(q.source(), tag) &&
      //	  !vqps.is_same_point(q.target(), tag) );
      if ( vqps.oriented_side(lq, tag) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    } else {
      if ( vpqr.oriented_side(lq, tag) != ON_POSITIVE_SIDE ) {
	lq = opposite_line(lq);
      }
    }

    // check if t is on the same side as the Voronoi vertices
    Oriented_side ot_lp = oriented_side_of_line(lp, t);
    Oriented_side ot_lq = oriented_side_of_line(lq, t);

    if ( ot_lp != ON_POSITIVE_SIDE || ot_lq != ON_POSITIVE_SIDE ) {
      return false;
    }

    Line lperp;

    Comparison_result res =
      compare_squared_distances_to_lines(t, lp, lq);

#if 0
    FT d2_t_from_lp = squared_distance_to_line(t.point(), lp);
    FT d2_t_from_lq = squared_distance_to_line(t.point(), lq);
      
    Comparison_result res =
      CGAL_NTS compare(d2_t_from_lp, d2_t_from_lq);
#endif

    if ( res == SMALLER ) {
      lperp = compute_perpendicular(lp, t);
    } else {
      lperp = compute_perpendicular(lq, t);
    }

    CGAL_precondition( ot_lp != ON_ORIENTED_BOUNDARY &&
		       ot_lq != ON_ORIENTED_BOUNDARY );

    // check of lperp separates the two Voronoi vertices
    Oriented_side opqr_perp = vpqr.oriented_side(lperp, tag);
    Oriented_side oqps_perp = vqps.oriented_side(lperp, tag);

    bool on_different_side =
      (opqr_perp == ON_POSITIVE_SIDE &&
       oqps_perp == ON_NEGATIVE_SIDE) ||
      (opqr_perp == ON_NEGATIVE_SIDE && 
       oqps_perp == ON_POSITIVE_SIDE);

    return ( on_different_side );
  }

  //------------------------------------------------------------------------



#if 0
  bool is_interior_in_conflict_none(const Site& p, const Site& q,
				    const Site& r, const Site& s,
				    const Site& t, Ring_tag tag) const
  {
    // THIS NEEDS TO BE CORRECTED
    return
      is_interior_in_conflict_none(p, q, r, s, t, Sqrt_field_tag());
  }



  bool
  is_interior_in_conflict_none(const Site& p, const Site& q,
			       const Site& r, const Site& s,
			       const Site& t, Sqrt_field_tag tag) const
  {
    // checks if interior of voronoi edge is in conflict if both extrema 
    // of the voronoi edge are not in conclift
    // return true if interior is in conflict; false otherwise

    if ( t.is_segment() ) { return false; }

    if ( p.is_point() && q.is_point() ) { return false; }

    if ( p.is_point() && q.is_segment() ) {
      if ( p.point() == q.source() || p.point() == q.target() ) {
	return false;
      }

      //      Line lq(q.segment());

      Line lq = compute_supporting_line(q.segment());
    
      Voronoi_vertex vpqr(p, q, r);
      Voronoi_vertex vqps(q, p, s);

      //      Line lperp = lq.perpendicular(t.point());
      Line lperp = compute_perpendicular(lq, t.point());

      //      Oriented_side op = lq.oriented_side(p.point());
      //      Oriented_side ot = lq.oriented_side(t.point());

      Oriented_side op = oriented_side_of_line(lq, p.point());
      Oriented_side ot = oriented_side_of_line(lq, t.point());

      bool on_same_side =
	((op == ON_POSITIVE_SIDE && ot == ON_POSITIVE_SIDE) ||
	 (op == ON_NEGATIVE_SIDE && ot == ON_NEGATIVE_SIDE));

#if 0
      FT d2_t_from_lq = squared_distance_to_line(t.point(), lq);
      FT d2_p_from_lq = squared_distance_to_line(p.point(), lq);

      Comparison_result res =
	CGAL_NTS compare(d2_t_from_lq, d2_p_from_lq);
#endif

      Comparison_result res =
	compare_squared_distances_to_line(lq, t.point(), p.point());

      //      Oriented_side opqr = lperp.oriented_side(vpqr.point());
      //      Oriented_side oqps = lperp.oriented_side(vqps.point());

      Oriented_side opqr = vpqr.oriented_side(lperp, tag);
      Oriented_side oqps = vqps.oriented_side(lperp, tag);

      bool on_different_side =
	((opqr == ON_POSITIVE_SIDE && oqps == ON_NEGATIVE_SIDE) ||
	 (opqr == ON_NEGATIVE_SIDE && oqps == ON_POSITIVE_SIDE));
    

      return ( on_same_side && (res == SMALLER) &&
	       on_different_side );
    } else if ( p.is_segment() && q.is_point() ) {
      return is_interior_in_conflict_none(q, p, s, r, t, tag);
    } else {
      // both p and q are segments

      Voronoi_vertex vpqr(p, q, r);
      Voronoi_vertex vqps(q, p, s);

      Line lp = compute_supporting_line(p.segment());
      Line lq = compute_supporting_line(q.segment());

      // first orient lp according to its Voronoi vertices
      if ( vpqr.is_same_point(p.source(), tag) ||
	   vpqr.is_same_point(p.target(), tag) ) {
	CGAL_assertion
	  ( !vqps.is_same_point(p.source(), tag) &&
	    !vqps.is_same_point(p.target(), tag) );
	//	if ( lp.oriented_side(vqps.point()) != ON_POSITIVE_SIDE ) {
	if ( vqps.oriented_side(lp, tag) != ON_POSITIVE_SIDE ) {
	  lp = opposite_line(lp);
	}
      } else {
	if ( vpqr.oriented_side(lp, tag) != ON_POSITIVE_SIDE ) {
	  lp = opposite_line(lp);
	}
      }

    // then orient lq according to its Voronoi vertices
      if ( vpqr.is_same_point(q.source(), tag) ||
	   vpqr.is_same_point(q.target(), tag) ) {
	CGAL_assertion
	  ( !vqps.is_same_point(q.source(), tag) &&
	    !vqps.is_same_point(q.target(), tag) );
	if ( vqps.oriented_side(lq, tag) != ON_POSITIVE_SIDE ) {
	  lq = opposite_line(lq);
	}
      } else {
	if ( vpqr.oriented_side(lq, tag) != ON_POSITIVE_SIDE ) {
	  lq = opposite_line(lq);
	}
      }

      // check if t is on the same side as the Voronoi vertices
      Oriented_side ot_lp = oriented_side_of_line(lp, t.point());
      Oriented_side ot_lq = oriented_side_of_line(lq, t.point());

      if ( ot_lp != ON_POSITIVE_SIDE ||
	   ot_lq != ON_POSITIVE_SIDE ) {
	return false;
      }

      Line lperp;

#if 0
      FT d2_t_from_lp = squared_distance_to_line(t.point(), lp);
      FT d2_t_from_lq = squared_distance_to_line(t.point(), lq);

      Comparison_result res =
	CGAL_NTS compare(d2_t_from_lp, d2_t_from_lq);
#endif

      Comparison_result res =
	compare_squared_distances_to_lines(t.point(), lp, lq);

      if ( res == SMALLER ) {
	lperp = compute_perpendicular(lp, t.point());
      } else {
	lperp = compute_perpendicular(lq, t.point());
      }

      CGAL_precondition( ot_lp != ON_ORIENTED_BOUNDARY &&
			 ot_lq != ON_ORIENTED_BOUNDARY );

      // check of lperp separates the two Voronoi vertices
      //      Oriented_side opqr_perp = lperp.oriented_side(vpqr.point());
      //      Oriented_side oqps_perp = lperp.oriented_side(vqps.point());

      Oriented_side opqr_perp = vpqr.oriented_side(lperp, tag);
      Oriented_side oqps_perp = vqps.oriented_side(lperp, tag);



      bool on_different_side =
	(opqr_perp == ON_POSITIVE_SIDE &&
	 oqps_perp == ON_NEGATIVE_SIDE) ||
	(opqr_perp == ON_NEGATIVE_SIDE && 
	 oqps_perp == ON_POSITIVE_SIDE);

      return ( on_different_side );
    }
  }
#endif

  //------------------------------------------------------------------------

public:
  inline
  bool operator()(const Site& p, const Site& q, const Site& r,
		  const Site& s, const Site& t, Sign sgn) const
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


  inline
  bool operator()(const Site& p, const Site& q, const Site& r,
		  const Site& t, Sign sgn) const
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

    bool p_is_endpoint = (p.point() == t.source() || p.point() == t.target());
    bool q_is_endpoint = (q.point() == t.source() || q.point() == t.target());

    return ( p_is_endpoint && q_is_endpoint );
  }

  inline
  bool operator()(const Site& p, const Site& q, const Site& t,
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

      //      Sign s1 = sign_of_determinant2x2(tp.x(), -tp.y(),
      //				       tq.y(), tq.x());

      CGAL_assertion( s1 != ZERO );
      return ( s1 == NEGATIVE );
    }

    bool bp = ( (p.point() == t.source()) || (p.point() == t.target()) );
    bool bq = ( (q.point() == t.source()) || (q.point() == t.target()) );

    return ( bp && bq );
  }

};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


template<class R, class Method_tag>
class Svd_infinite_edge_interior_2
{
public:
  typedef typename R::Point_2     Point;
  typedef typename R::Segment_2   Segment;
  typedef typename R::Line_2      Line;
  typedef typename R::Vector_2    Vector;
  typedef typename CGAL::Svd_voronoi_vertex_2<R>  Voronoi_vertex;
  typedef typename R::Site_2      Site;
  typedef typename R::FT                                      FT;

public:
  inline
  bool operator()(const Site& q, const Site& s, const Site& r,
		  const Site& t, Sign sgn)
  {
    if ( t.is_segment() ) {
#if PRED_PRINT
      std::cout << false << std::endl;
      return false;
#endif
    }

#if 0
    if ( q.is_segment() ) {
      // in this case r and s must be endpoints of q
      return ( sgn == NEGATIVE );
    }
#endif

    return ( sgn == NEGATIVE );
  }

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


template<class R, class Method_tag>
class Svd_is_degenerate_edge_test_2
{
public:
  typedef typename R::Point_2     Point;
  typedef typename R::Segment_2   Segment;
  typedef typename R::Line_2      Line;
  typedef typename R::Vector_2    Vector;
  typedef typename CGAL::Svd_voronoi_vertex_2<R>  Voronoi_vertex;
  typedef typename R::Site_2      Site;
  typedef typename R::FT          FT;

public:
  inline
  bool operator()(const Site& p, const Site& q, const Site& r,
		  const Site& s)
  {
    return false;
  }

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H

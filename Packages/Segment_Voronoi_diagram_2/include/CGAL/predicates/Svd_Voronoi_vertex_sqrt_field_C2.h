// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/predicates/Svd_Voronoi_vertex_sqrt_field_C2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_SQRT_FIELD_C2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VERTEX_SQRT_FIELD_C2_H



#include <CGAL/enum.h>
#include <CGAL/predicates/Svd_basic_predicates_C2.h>



CGAL_BEGIN_NAMESPACE




template<class K>
class Svd_voronoi_vertex_sqrt_field_C2
  : public Svd_basic_predicates_C2<K>
{
public:
  typedef Svd_basic_predicates_C2<K> Base;

  typedef enum {PPP = 0, PPS, PSS, SSS} vertex_t;
  struct PPP_Type {};
  struct PPS_Type {};
  struct PSS_Type {};
  struct SSS_Type {};

  typedef typename Base::Point_2             Point_2;
  typedef typename Base::Segment_2           Segment_2;
  typedef typename Base::Line_2              Line_2;
  typedef typename Base::Site_2              Site_2;
  typedef typename Base::FT                  FT;
  typedef typename Base::RT                  RT;

  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

private:
  //--------------------------------------------------------------------------

  void
  compute_ppp(const Point_2& p, const Point_2& q, const Point_2& r)
  {
    v_type = PPP;

    FT np = CGAL::square(p.x()) + CGAL::square(p.y());
    FT nq = CGAL::square(q.x()) + CGAL::square(q.y());
    FT nr = CGAL::square(r.x()) + CGAL::square(r.y());

    ux = np * (q.y() - r.y()) + nq * (r.y() - p.y()) + nr * (p.y() - q.y());
    uy = -(np * (q.x() - r.x()) + nq * (r.x() - p.x()) + nr * (p.x() - q.x()));
    uz = FT(2) * ( (q.x() * r.y() - r.x() * q.y()) +
		   (r.x() * p.y() - p.x() * r.y()) +
		   (p.x() * q.y() - q.x() * p.y()) );
  }

  //--------------------------------------------------------------------------

  void
  compute_pss(const Point_2& p, const Segment_2& q, const Segment_2& r)
  {
    v_type = PSS;

    bool pq = (p == q.source() || p == q.target());
    bool pr = (p == r.source() || p == r.target());

    if ( pq && pr ) {
      ux = p.x();
      uy = p.y();
      uz = FT(1);
      return;
    }


    FT a1, b1, c1, a2, b2, c2;
    compute_supporting_line(q, a1, b1, c1);
    compute_supporting_line(r, a2, b2, c2);

    FT c1_ = a1 * p.x() + b1 * p.y() + c1;
    FT c2_ = a2 * p.x() + b2 * p.y() + c2;

    if ( pq ) {
      c1_ = FT(0);
    }

    if ( pr ) {
      c2_ = FT(0);
    }

    Sign sgn_c1_ = CGAL::sign(c1_);
    Sign sgn_c2_ = CGAL::sign(c2_);

    if ( sgn_c1_ == NEGATIVE ) {
      a1 = -a1;  b1 = -b1;  c1_ = -c1_;
    } else if ( sgn_c1_ == ZERO ) {

      CGAL_assertion( pq );
      if ( p == q.target() ) {
	a1 = -a1;  b1 = -b1;  c1_ = -c1_;
      }
    }

    if ( sgn_c2_ == NEGATIVE ) {
      a2 = -a2;  b2 = -b2;  c2_ = -c2_;
    } else if ( sgn_c2_ == ZERO ) {

      CGAL_assertion( pr );
      if ( p == r.source() ) {
	a2 = -a2;  b2 = -b2;  c2_ = -c2_;
      }
    }

    if ( pq ) {
      FT J = a1 * c2_;
      FT I = b1 * c2_;

      FT n1 = CGAL::square(a1) + CGAL::square(b1);
      FT n2 = CGAL::square(a2) + CGAL::square(b2);

      FT D1D2 = n1 * n2;

      uz = -a1 * a2 - b1 * b2 + CGAL::sqrt(D1D2);

      ux = J + p.x() * uz;
      uy = I + p.y() * uz;

    } else if ( pr ) {
      FT J = a2 * c1_;
      FT I = b2 * c1_;

      FT n1 = CGAL::square(a1) + CGAL::square(b1);
      FT n2 = CGAL::square(a2) + CGAL::square(b2);

      FT D1D2 = n1 * n2;

      uz = -a1 * a2 - b1 * b2 + CGAL::sqrt(D1D2);

      ux = J + p.x() * uz;
      uy = I + p.y() * uz;

    } else {
      Line_2 lq(a1, b1, c1_);
      Line_2 lr(a2, b2, c2_);
      compute_pll(p, lq, lr);
    }
  }


  void
  compute_pll(const Point_2& p, const Line_2& lq, const Line_2& lr)
  {
    FT a1 = lq.a(), b1 = lq.b(), c1_ = lq.c();
    FT a2 = lr.a(), b2 = lr.b(), c2_ = lr.c();

    CGAL_precondition( c1_ >= FT(0) );
    CGAL_precondition( c2_ >= FT(0) );

    FT n1 = CGAL::square(a1) + CGAL::square(b1);
    FT n2 = CGAL::square(a2) + CGAL::square(b2);

    FT I = b1 * c2_ + b2 * c1_;
    FT J = a1 * c2_ + a2 * c1_;

    FT c1c2 = FT(2) * c1_ * c2_;
    FT a1a2 = a1 * a2;
    FT b1b2 = b1 * b2;

    FT D1D2 = n1 * n2;

    // compute sigma
    FT sigma_expr = b1 * CGAL::sqrt(n2) - b2 * CGAL::sqrt(n1);
    Sign s_sigma = CGAL::sign(sigma_expr);

    int i_sigma(s_sigma);
    FT sigma(i_sigma);

    // compute rho
    FT rho_expr = a1 * CGAL::sqrt(n2) - a2 * CGAL::sqrt(n1);
    Sign s_rho = CGAL::sign(rho_expr);

    int i_rho(s_rho);
    FT rho(i_rho);

    FT sqrt_D1D2 = CGAL::sqrt(D1D2);

    FT A = a1a2 - b1b2;
    FT u1 = c1c2 * (sqrt_D1D2 + A);
    FT u2 = c1c2 * (sqrt_D1D2 - A);

    uz = -a1a2 - b1b2 + CGAL::sqrt(D1D2);

    ux = J + p.x() * uz + sigma * CGAL::sqrt(u1);
    uy = I + p.y() * uz - rho * CGAL::sqrt(u2);
  }


  //--------------------------------------------------------------------------


  void
  compute_pps(const Point_2& p, const Point_2& q, const Segment_2& r)
  {
    v_type = PPS;

    FT a, b, c;
    compute_supporting_line(r, a, b, c);

    FT c_ = a * p.x() + b * p.y() + c;
    FT cq_ = a * q.x() + b * q.y() + c;


    if ( p == r.source() || p == r.target() ) {
      c_ = FT(0);
    }
    if ( q == r.source() || q == r.target() ) {
      cq_ = FT(0);
    }

    Sign s = CGAL::sign(c_);

    if ( s == NEGATIVE ) {
      a = -a;  b = -b;  c = -c;  c_ = -c_;  cq_ = -cq_;
    } else if ( s == ZERO ) {
      Sign s1 = CGAL::sign(cq_);

      CGAL_assertion( s1 != ZERO );
      if ( s1 == NEGATIVE ) {
	a = -a;  b = -b;  c = -c;  c_ = -c_;  cq_ = -cq_;
      }
    }

    FT nl = CGAL::square(a) + CGAL::square(b);

    FT x_ = q.x() - p.x();
    FT y_ = q.y() - p.y();
    FT n_ = CGAL::square(x_) + CGAL::square(y_);


    Comparison_result res = CGAL::compare( c_, cq_ );

    if ( res == EQUAL ) {
      FT e1 = CGAL::square(c_);
      FT J = nl * (a * n_ + FT(4) * c_ * x_) - FT(4) * a * e1;
      FT I = nl * (b * n_ + FT(4) * c_ * y_) - FT(4) * b * e1;
      FT X = FT(8) * nl * c_;

      ux = J + p.x() * X;
      uy = I + p.y() * X;
      uz = X;
      return;
    }


    FT e1 = a * x_ + b * y_;
    FT e2 = b * x_ - a * y_;
    FT e3 = n_ * e1;
    FT e4 = FT(2) * c_ * e2;

    FT X = FT(2) * CGAL::square(e1);
    FT I = b * e3 + x_ * e4;
    FT J = a * e3 - y_ * e4;
    FT sqrt_S = CGAL::sqrt(n_ * nl * c_ * cq_);

    ux = J + p.x() * X - FT(2) * y_ * sqrt_S;
    uy = I + p.y() * X + FT(2) * x_ * sqrt_S;
    uz = X;
  }


  //--------------------------------------------------------------------------

  bool
  is_consistent(const Segment_2& p, const Segment_2& q,
		const Segment_2& r, FT a[], FT b[], FT c[]) const
  {
    Line_2 l[3] = {Line_2(a[0], b[0], c[0]), Line_2(a[1], b[1], c[1]), 
		   Line_2(a[2], b[2], c[2])};

    int num_oriented(0);

    if ( is_on_positive_halfspace(l[0], q) ||
	 is_on_positive_halfspace(l[0], r) ) {
      num_oriented++;
    }

    if ( is_on_positive_halfspace(l[1], p) ||
	 is_on_positive_halfspace(l[1], r) ) {
      num_oriented++;
    }

    if ( is_on_positive_halfspace(l[2], p) ||
	 is_on_positive_halfspace(l[2], q) ) {
      num_oriented++;
    }

    return ( num_oriented >= 2 );
  }

  void
  orient_lines(const Segment_2& p, const Segment_2& q,
	       const Segment_2& r, FT a[], FT b[], FT c[]) const 
  {
    Line_2 l[3];
    l[0] = compute_supporting_line(p);
    l[1] = compute_supporting_line(q);
    l[2] = compute_supporting_line(r);
    
    bool is_oriented[3] = {false, false, false};

    if ( is_on_positive_halfspace(l[0], q) ||
	 is_on_positive_halfspace(l[0], r) ) {
      is_oriented[0] = true;
    } else {
      l[0] = opposite_line(l[0]);
      if ( is_on_positive_halfspace(l[0], q) ||
	   is_on_positive_halfspace(l[0], r) ) {
	is_oriented[0] = true;
      } else {
	l[0] = opposite_line(l[0]);
      }
    }

    if ( is_on_positive_halfspace(l[1], p) ||
	 is_on_positive_halfspace(l[1], r) ) {
      is_oriented[1] = true;
    } else {
      l[1] = opposite_line(l[1]);
      if ( is_on_positive_halfspace(l[1], p) ||
	   is_on_positive_halfspace(l[1], r) ) {
	is_oriented[1] = true;
      } else {
	l[1] = opposite_line(l[1]);
      }
    }

    if ( is_on_positive_halfspace(l[2], p) ||
	 is_on_positive_halfspace(l[2], q) ) {
      is_oriented[2] = true;
    } else {
      l[2] = opposite_line(l[2]);
      if ( is_on_positive_halfspace(l[2], p) ||
	   is_on_positive_halfspace(l[2], q) ) {
	is_oriented[2] = true;
      } else {
	l[2] = opposite_line(l[2]);
      }
    }

    if ( is_oriented[0] && is_oriented[1] && is_oriented[2] ) {
      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
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

    FT sqrt_d[3];
    for (int i = 0; i < 3; i++) {
      FT d1 = CGAL::square(l[i].a()) + CGAL::square(l[i].b());
      sqrt_d[i] = CGAL::sqrt(d1);
    }

    FT z[3];
    for (int i = 0; i < 3; i++) {
      z[i] = l[(i+1)%3].a() * l[(i+2)%3].b()
	- l[(i+2)%3].a() * l[(i+1)%3].b();
    }


    FT vz = z[0] * sqrt_d[0] + z[1] * sqrt_d[1] + z[2] * sqrt_d[2];

    Sign s_minus_vz = CGAL::sign(vz);

    CGAL_assertion( s_minus_vz != ZERO );

    if ( s_minus_vz == NEGATIVE ) {
      l[i_no] = opposite_line(l[i_no]);

      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
      return;
    }

    // now we have to check if the other orientation of l[i_no]
    // corresponds to a CCW triangle as well.
    z[(i_no+1)%3] = -z[(i_no+1)%3];
    z[(i_no+2)%3] = -z[(i_no+2)%3];

    vz = z[0] * sqrt_d[0] + z[1] * sqrt_d[1] + z[2] * sqrt_d[2];

    Sign s_minus_vz_2 = CGAL::sign(vz);

    CGAL_assertion( s_minus_vz_2 != ZERO );

    if ( s_minus_vz_2 == NEGATIVE ) {
      // the other orientation does not correspond to a CCW triangle.
      for (int i = 0; i < 3; i++) {
	a[i] = l[i].a();
	b[i] = l[i].b();
	c[i] = l[i].c();
      }
      return;
    }

    // now compute the Voronoi vertex;
    for (int i = 0; i < 3; i++) {
      a[i] = l[i].a();
      b[i] = l[i].b();
      c[i] = l[i].c();
    }

    FT x[3], y[3], w[3];
    for (int i = 0; i < 3; i++) {
      x[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      y[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      w[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);
    }

    FT vx = FT(0), vy = FT(0), vw = FT(0);

    for (int i = 0; i < 3; i++) {
      vx += x[i] * sqrt_d[i];
      vy += y[i] * sqrt_d[i];
      vw += w[i] * sqrt_d[i];
    }

    Sign s_vw = CGAL::sign(vw);

    FT dist =
      a[(i_no+1)%3] * vx + b[(i_no+1)%3] * vy + c[(i_no+1)%3] * vw;
    

    Sign sgn_dist = Sign(s_vw * CGAL::sign(dist));

    CGAL_assertion( sgn_dist != ZERO );

    if ( sgn_dist == NEGATIVE ) {
      a[i_no] = -a[i_no];
      b[i_no] = -b[i_no];
      c[i_no] = -c[i_no];
    }
  }



  void
  compute_sss(const Segment_2& p, const Segment_2& q, const Segment_2& r)
  {
    v_type = SSS;

    FT a[3], b[3], c[3];
    FT cx[3], cy[3], cz[3], sqrt_D[3];

    orient_lines(p, q, r, a, b, c);

    for (int i = 0; i < 3; i++) {
      cx[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      cy[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      cz[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);

      FT d = CGAL::square(a[i]) + CGAL::square(b[i]);
      sqrt_D[i] = CGAL::sqrt(d);
    }

    ux = cx[0] * sqrt_D[0] + cx[1] * sqrt_D[1] + cx[2] * sqrt_D[2];
    uy = cy[0] * sqrt_D[0] + cy[1] * sqrt_D[1] + cy[2] * sqrt_D[2];
    uz = cz[0] * sqrt_D[0] + cz[1] * sqrt_D[1] + cz[2] * sqrt_D[2];
  }

  //--------------------------------------------------------------------------

  void
  compute_vertex(const Site_2& s1, const Site_2& s2, const Site_2& s3)
  {
    if ( s1.is_point() && s2.is_point() && s3.is_point() ) {
      compute_ppp(s1.point(), s2.point(), s3.point());

    } else if ( s1.is_segment() && s2.is_point() && s3.is_point() ) {
      compute_vertex(s2, s3, s1);
      pps_idx = 1;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
      pps_idx = 2;

    } else if ( s1.is_point() && s2.is_point() && s3.is_segment() ) {
      compute_pps(s1.point(), s2.point(), s3.segment());
      pps_idx = 0;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_segment() ) {
      compute_pss(s1.point(), s2.segment(), s3.segment());
    } else if ( s1.is_segment() && s2.is_point() && s3.is_segment() ) {
      compute_vertex(s2, s3, s1);
    } else if ( s1.is_segment() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
    } else {
      compute_sss(s1.segment(), s2.segment(), s3.segment());
    }

  }


  //--------------------------------------------------------------------------

  bool are_identical(const Point_2& p, const Point_2& q) const
  {
    return (p == q);
  }

  //--------------------------------------------------------------------------

  bool is_endpoint_of(const Point_2& p, const Segment_2& s) const
  {
    return ( are_identical(p, s.source()) ||
	     are_identical(p, s.target()) );
  }
  

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                              the incircle test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a point
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------

  Sign
  check_easy_degeneracies(const Point_2& t, PPS_Type,
			  bool& use_result) const
  {
    use_result = false;
    if (  ( p_.is_point() && are_identical(p_.point(),t) ) ||
	  ( q_.is_point() && are_identical(q_.point(),t) ) ||
	  ( r_.is_point() && are_identical(r_.point(),t) )  ) {
      use_result = true;
      return ZERO;
    }

    if (  ( p_.is_segment() && is_endpoint_of(t, p_.segment()) ) || 
	  ( q_.is_segment() && is_endpoint_of(t, q_.segment()) ) ||
	  ( r_.is_segment() && is_endpoint_of(t, r_.segment()) )  ) {
      use_result = true;
      return POSITIVE;
    }

    return ZERO;
  }

  Sign
  check_easy_degeneracies(const Point_2& t, PSS_Type,
			  bool& use_result) const
  {
    return check_easy_degeneracies(t, PPS_Type(), use_result);
  }

  Sign
  check_easy_degeneracies(const Point_2& t, SSS_Type,
			  bool& use_result) const
  {
    use_result = false;

    // ADD THE CASES WHERE t IS AN ENDPOINT OF ONE OF THE SEGMENTS
    return ZERO;
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Point_2& t, PPP_Type) const
  {
    Oriented_side os =
      side_of_oriented_circle(p_.point(), q_.point(), r_.point(), t);
    if ( os == ON_POSITIVE_SIDE ) { return NEGATIVE; }
    if ( os == ON_NEGATIVE_SIDE ) { return POSITIVE; }
    return ZERO;
  }

  //--------------------------------------------------------------------------

  template<class Type>
  Sign incircle(const Point_2& t, Type type) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    FT r2 = squared_radius();

    FT d2 = CGAL::square(x() - t.x()) +
      CGAL::square(y() - t.y());

    return Sign( CGAL::compare(d2, r2) );
  }


  //--------------------------------------------------------------------------

  Sign incircle(const Point_2& t) const 
  {
    if ( is_degenerate_Voronoi_circle() ) {
      return POSITIVE;
    }

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle(t, PPP_Type());
      break;
    case PPS:
      s = incircle(t, PPS_Type());
      break;
    case PSS:
      s = incircle(t, PSS_Type());
      break;
    case SSS:
      s = incircle(t, SSS_Type());
      break;
    }

    return s;
  }

  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a segment
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------


  Oriented_side
  oriented_side(const Line_2& l, const Point_2& p) const
  {
    Line_2 l1(l.b(), -l.a(), l.a() * y() - l.b() * x());

    return oriented_side_of_line(l1, p);
  }


  Sign incircle(const Line_2& l) const
  {
    FT r2 = squared_radius();

    FT n2 = CGAL::square(l.a()) + CGAL::square(l.b());

    FT d2 = CGAL::square(l.a() * x() + l.b() * y() + l.c());

    return Sign( CGAL::compare(d2, r2 * n2) );
  }



  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  Sign incircle1(const Segment_2& t) const
  {
    if ( v_type == PPP || v_type == PPS ) {
      if (  p_.is_point() && q_.is_point() &&
	    is_endpoint_of(p_.point(), t) &&
	    is_endpoint_of(q_.point(), t)  ) {
	return NEGATIVE;
      }

      if (  p_.is_point() && r_.is_point() &&
	    is_endpoint_of(p_.point(), t) &&
	    is_endpoint_of(r_.point(), t)  ){
	return NEGATIVE;
      }

      if (  q_.is_point() && r_.is_point() &&
	    is_endpoint_of(q_.point(), t) &&
	    is_endpoint_of(r_.point(), t)  ){
	return NEGATIVE;
      }
    }

    Sign d1, d2;
    if (  ( p_.is_point() && are_identical(p_.point(), t.source()) ) ||
	  ( q_.is_point() && are_identical(q_.point(), t.source()) ) ||
	  ( r_.is_point() && are_identical(r_.point(), t.source()) )  ) {
      d1 = ZERO;
    } else {
      d1 = incircle(t.source());
    }
    if ( d1 == NEGATIVE ) { return NEGATIVE; }

    if (  ( p_.is_point() && are_identical(p_.point(), t.target()) ) ||
	  ( q_.is_point() && are_identical(q_.point(), t.target()) ) ||
	  ( r_.is_point() && are_identical(r_.point(), t.target()) )  ) {
      d2 = ZERO;
    } else {
      d2 = incircle(t.target());
    }
    if ( d2 == NEGATIVE ) { return NEGATIVE; }


    Line_2 l = compute_supporting_line(t);
    Sign sl = incircle(l);

    // this is the old code
    if ( sl == POSITIVE ) { return sl; }

    Oriented_side os1 = oriented_side(l, t.source());
    Oriented_side os2 = oriented_side(l, t.target());

    if ( sl == ZERO ) {
      if ( (os1 == ON_POSITIVE_SIDE && os2 == ON_NEGATIVE_SIDE) ||
	   (os1 == ON_NEGATIVE_SIDE && os2 == ON_POSITIVE_SIDE) ) {
	return ZERO;
      }
      return POSITIVE;
    }

    return (os1 == os2) ? POSITIVE : NEGATIVE;
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Segment_2& t) const 
  {
    if ( is_degenerate_Voronoi_circle() ) {
      // case 1: the new segment is not adjacent to the center of the
      //         degenerate Voronoi circle
      if (  !are_identical( p_ref(), t.source() ) &&
	    !are_identical( p_ref(), t.target() )  ) {
	return POSITIVE;
      }

      if (  ( r_.is_point() && are_identical(p_ref(),r_.point()) ) ||
	    ( q_.is_point() && are_identical(p_ref(),q_.point()) ) ||
	    ( p_.is_point() && are_identical(p_ref(),p_.point()) )  ) {
	Point_2 pr;
	Segment_2 sp, sq;
	if ( p_.is_point() ) {
	  CGAL_assertion( q_.is_segment() && r_.is_segment() );
	  pr = p_.point();
	  sp = q_.segment();
	  sq = r_.segment();
	} else if ( q_.is_point() ) {
	  CGAL_assertion( r_.is_segment() && p_.is_segment() );
	  pr = q_.point();
	  sp = r_.segment();
	  sq = p_.segment();
	} else {
	  CGAL_assertion( p_.is_segment() && q_.is_segment() );
	  pr = r_.point();
	  sp = p_.segment();
	  sq = q_.segment();
	}

	Point_2 pq = sq.source(), pp = sp.source(), pt = t.source();

	if ( are_identical(sp.source(),pr) ) { pp = sp.target(); }
	if ( are_identical(sq.source(),pr) ) { pq = sq.target(); }
	if ( are_identical( t.source(),pr) ) { pt =  t.target(); }

	if ( CGAL::orientation(pr, pp, pt) == LEFT_TURN &&
	     CGAL::orientation(pr, pq, pt) == RIGHT_TURN ) {
	  return NEGATIVE;
	}
	return ZERO;
      }
    } // if ( is_degenerate_Voronoi_circle() )

    Sign s = incircle1(t);

    return s;
  }

  //--------------------------------------------------------------------------
  //  subpredicates for the incircle test
  //--------------------------------------------------------------------------


public:
  bool is_degenerate_Voronoi_circle() const
  {
    if ( v_type != PSS ) { return false; }

    if ( p_.is_point() ) {
      return ( is_endpoint_of(p_.point(), q_.segment()) &&
	       is_endpoint_of(p_.point(), r_.segment()) );
    } else if ( q_.is_point() ) {
      return ( is_endpoint_of(q_.point(), p_.segment()) &&
	       is_endpoint_of(q_.point(), r_.segment()) );
    } else {
      CGAL_assertion( r_.is_point() );
      return ( is_endpoint_of(r_.point(), p_.segment()) &&
	       is_endpoint_of(r_.point(), q_.segment()) );
    }
  }


  //--------------------------------------------------------------------------


public:
  // THIS METHOD SHOULD BE PRIVATE AND SHOULD NOT BE USED IN THE
  // PREDICATES (MAYBE I AM WRONG ON THIS????) IN ANY CASE CALLS TO
  // THIS METHOD SHOULD BE ELIMINATED AS MUCH AS POSSIBLE; ALSO I NEED
  // TO WRITE A SPECIAL VERSION FOR THE RING_TAG CASE
#if 0
  bool is_same_point(const Point_2& p) const
  {
    Comparison_result res = CGAL::compare(x(), p.x());
    if ( res != EQUAL ) { return false; }

    return (CGAL::compare(y(), p.y()) == EQUAL);
  }
#endif

private:

  //--------------------------------------------------------------------------
  //  the reference point (valid if v_type != SSS)
  //--------------------------------------------------------------------------

  const Point_2& p_ref() const
  {
    CGAL_precondition ( v_type != SSS );

    if ( v_type == PPS ) {
      if ( pps_idx == 0 ) { return p_.point(); }
      if ( pps_idx == 1 ) { return q_.point(); }
      return r_.point();
    }

    if ( p_.is_point() ) {
      return p_.point();
    } else if ( q_.is_point() ) {
      return q_.point();
    } else {
      CGAL_assertion( r_.is_point() );
      return r_.point();
    }
  }


public:
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                           access methods
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  FT x() const { return hx() / hw(); }
  FT y() const { return hy() / hw(); }

  FT hx() const {
    return ux;
  }

  FT hy() const {
    return uy;
  }

  FT hw() const {
    return uz;
  }

  FT squared_radius() const {
    switch (v_type) {
    case PPP:    case PPS:    case PSS:
      {
	FT dx2 = CGAL::square(x() - p_ref().x());
	FT dy2 = CGAL::square(y() - p_ref().y());
	return dx2 + dy2;
      }
      break;
    case SSS:
      {
	Line_2 l = compute_supporting_line(p_.segment());
	Homogeneous_point_2 q = compute_projection(l, point());

	FT dx2 = CGAL::square(x() - q.x());
	FT dy2 = CGAL::square(y() - q.y());
	return dx2 + dy2;
      }
      break;
    default:
      return FT(0);
    }
  }


  Point_2 point() const {
    if ( is_degenerate_Voronoi_circle() ) {
      return p_ref();
    }
    
    return Point_2(x(), y());
  }


  Point_2 degenerate_point() const
  {
    CGAL_precondition( is_degenerate_Voronoi_circle() );
    return p_ref();
  }


  typename K::Circle_2 circle() const
  {
    typedef typename K::Circle_2  K_Circle_2;
    return K_Circle_2(point(), squared_radius());
  }

  vertex_t type() const { return v_type; }

public:
  Svd_voronoi_vertex_sqrt_field_C2(const Site_2& p,
				   const Site_2& q,
				   const Site_2& r)
    : p_(p), q_(q), r_(r)
  {
    compute_vertex(p, q, r);
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Site_2& t) const 
  {
    if ( t.is_point() ) {
      return incircle(t.point());
    }
    return incircle(t.segment());
  }

  //--------------------------------------------------------------------------


  Orientation orientation(const Line_2& l) const 
  {
    Sign s = CGAL::sign(l.a() * x() + l.b() * y() + l.c());

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }

  Oriented_side oriented_side(const Line_2& l) const
  {
    Orientation o = orientation(l);

    if ( o == COLLINEAR ) { return ON_ORIENTED_BOUNDARY; }
    return ( o == LEFT_TURN ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  }

  //--------------------------------------------------------------------------

private:
  const Site_2& p_, q_, r_;

  vertex_t v_type;

  // index that indicates the refence point for the case PPS
  short pps_idx;

  FT ux, uy, uz;
};




CGAL_END_NAMESPACE



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_VEFTEX_SQRT_FIELD_C2_H

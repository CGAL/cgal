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
// file          : include/CGAL/predicates/Svd_Voronoi_vertex_ring_C2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_RING_C2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_RING_C2_H


#include <CGAL/enum.h>
#include <CGAL/predicates/Svd_basic_predicates_C2.h>



CGAL_BEGIN_NAMESPACE



template<class K>
class Svd_voronoi_vertex_ring_C2
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

  typedef typename Base::Sqrt_1              Sqrt_1;
  typedef typename Base::Sqrt_2              Sqrt_2;
  typedef typename Base::Sqrt_3              Sqrt_3;

  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

private:
  //--------------------------------------------------------------------------

  void
  compute_ppp(const Point_2& p, const Point_2& q, const Point_2& r)
  {
    v_type = PPP;

    RT np = CGAL::square(p.x()) + CGAL::square(p.y());
    RT nq = CGAL::square(q.x()) + CGAL::square(q.y());
    RT nr = CGAL::square(r.x()) + CGAL::square(r.y());

    ux_ppp = 
      np * (q.y() - r.y()) + nq * (r.y() - p.y()) + nr * (p.y() - q.y());
    uy_ppp =
      -(np * (q.x() - r.x()) + nq * (r.x() - p.x()) + nr * (p.x() - q.x()));
    uz_ppp = RT(2) * ( (q.x() * r.y() - r.x() * q.y()) +
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
      Sqrt_1 One(RT(1), RT(0), RT(0));

      ux = Sqrt_3(p.x() * One);
      uy = Sqrt_3(p.y() * One);
      uz = Sqrt_3(One);
      return;
    }


    RT a1, b1, c1, a2, b2, c2;
    compute_supporting_line(q, a1, b1, c1);
    compute_supporting_line(r, a2, b2, c2);

    RT c1_ = a1 * p.x() + b1 * p.y() + c1;
    RT c2_ = a2 * p.x() + b2 * p.y() + c2;

    if ( pq ) {
      c1_ = RT(0);
    }

    if ( pr ) {
      c2_ = RT(0);
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
      RT J = a1 * c2_;
      RT I = b1 * c2_;

      RT n1 = CGAL::square(a1) + CGAL::square(b1);
      RT n2 = CGAL::square(a2) + CGAL::square(b2);

      RT D1D2 = n1 * n2;

      Sqrt_1 Zero(RT(0), RT(0), D1D2);

      Sqrt_1 vz(-a1 * a2 - b1 * b2, RT(1), D1D2);

      ux = Sqrt_3(J + p.x() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uy = Sqrt_3(I + p.y() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uz = Sqrt_3(vz, Zero, Zero, Zero, Zero, Zero);
    } else if ( pr ) {
      RT J = a2 * c1_;
      RT I = b2 * c1_;

      RT n1 = CGAL::square(a1) + CGAL::square(b1);
      RT n2 = CGAL::square(a2) + CGAL::square(b2);

      RT D1D2 = n1 * n2;

      Sqrt_1 Zero(RT(0), RT(0), D1D2);

      Sqrt_1 vz(-a1 * a2 - b1 * b2, RT(1), D1D2);

      ux = Sqrt_3(J + p.x() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uy = Sqrt_3(I + p.y() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uz = Sqrt_3(vz, Zero, Zero, Zero, Zero, Zero);
    } else {
      Line_2 lq(a1, b1, c1_);
      Line_2 lr(a2, b2, c2_);
      compute_pll(p, lq, lr);
    }
  }


  void
  compute_pll(const Point_2& p, const Line_2& lq, const Line_2& lr)
  {
    RT a1 = lq.a(), b1 = lq.b(), c1_ = lq.c();
    RT a2 = lr.a(), b2 = lr.b(), c2_ = lr.c();

    CGAL_precondition( c1_ >= RT(0) );
    CGAL_precondition( c2_ >= RT(0) );

    RT n1 = CGAL::square(a1) + CGAL::square(b1);
    RT n2 = CGAL::square(a2) + CGAL::square(b2);

    RT I = b1 * c2_ + b2 * c1_;
    RT J = a1 * c2_ + a2 * c1_;
    
    RT c1c2 = RT(2) * c1_ * c2_;
    RT a1a2 = a1 * a2;
    RT b1b2 = b1 * b2;

    RT D1D2 = n1 * n2;

    Sqrt_1 Zero(RT(0), RT(0), D1D2);
    Sqrt_1 One(RT(1), RT(0), D1D2);


    Sign s1, s2;

    // compute sigma
    Sign s_sigma(ZERO);
    s1 = CGAL::sign(b1);
    s2 = CGAL::sign(-b2);
    if ( s1 == ZERO ) {
      s_sigma = s2;
    } else if ( s2 == ZERO ) {
      s_sigma = s1;
    } else if ( s1 == s2 ) {
      s_sigma = s1;
    } else {
      RT e = CGAL::square(b1) * n2 - CGAL::square(b2) * n1;
      s_sigma = Sign(s1 * CGAL::sign(e));
    }

    Sqrt_1 sigma = Zero;
    if ( s_sigma == POSITIVE ) {
      sigma = One;
    } else if ( s_sigma == NEGATIVE ) {
      sigma = -One;
    }

    // compute rho
    Sign s_rho(ZERO);
    s1 = CGAL::sign(a1);
    s2 = CGAL::sign(-a2);
    if ( s1 == ZERO ) {
      s_rho = s2;
    } else if ( s2 == ZERO ) {
      s_rho = s1;
    } else if ( s1 == s2 ) {
      s_rho = s1;
    } else {
      RT e = CGAL::square(a1) * n2 - CGAL::square(a2) * n1;
      s_rho = Sign(s1 * CGAL::sign(e));
    }

    
    Sqrt_1 rho = Zero;
    if ( s_rho == POSITIVE ) {
      rho = One;
    } else if ( s_rho == NEGATIVE ) {
      rho = -One;
    }

    
    Sqrt_1 vz(-a1a2 - b1b2, RT(1), D1D2);

    RT A = a1a2 - b1b2;
    Sqrt_1 u1( c1c2 * A, c1c2, D1D2);
    Sqrt_1 u2(-c1c2 * A, c1c2, D1D2);


    ux = Sqrt_3(J + p.x() * vz, sigma, Zero, Zero, u1, u2);
    uy = Sqrt_3(I + p.y() * vz, Zero, -rho, Zero, u1, u2);
    uz = Sqrt_3(vz, Zero, Zero, Zero, u1, u2);
  }


  //--------------------------------------------------------------------------


  void
  compute_pps(const Point_2& p, const Point_2& q, const Segment_2& r)
  {
    v_type = PPS;

    RT a, b, c;
    compute_supporting_line(r, a, b, c);

    RT c_ = a * p.x() + b * p.y() + c;
    RT cq_ = a * q.x() + b * q.y() + c;


    if ( p == r.source() || p == r.target() ) {
      c_ = RT(0);
    }
    if ( q == r.source() || q == r.target() ) {
      cq_ = RT(0);
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

    RT nl = CGAL::square(a) + CGAL::square(b);

    RT x_ = q.x() - p.x();
    RT y_ = q.y() - p.y();
    RT n_ = CGAL::square(x_) + CGAL::square(y_);


    Comparison_result res = CGAL::compare( c_, cq_ );

    if ( res == EQUAL ) {
      RT e1 = CGAL::square(c_);
      RT J = nl * (a * n_ + RT(4) * c_ * x_) - RT(4) * a * e1;
      RT I = nl * (b * n_ + RT(4) * c_ * y_) - RT(4) * b * e1;
      RT X = RT(8) * nl * c_;

      ux_pps = Sqrt_1(J + p.x() * X);
      uy_pps = Sqrt_1(I + p.y() * X);
      uz_pps = Sqrt_1(X);
      return;
    }


    RT e1 = a * x_ + b * y_;
    RT e2 = b * x_ - a * y_;
    RT e3 = n_ * e1;
    RT e4 = RT(2) * c_ * e2;

    RT X = RT(2) * CGAL::square(e1);
    RT I = b * e3 + x_ * e4;
    RT J = a * e3 - y_ * e4;
    RT S = n_ * nl * c_ * cq_;

    ux_pps = Sqrt_1(J + p.x() * X, RT(-2) * y_, S);
    uy_pps = Sqrt_1(I + p.y() * X, RT( 2) * x_, S);
    uz_pps = Sqrt_1(X, RT(0), S);
  }


  //--------------------------------------------------------------------------

  bool
  is_consistent(const Segment_2& p, const Segment_2& q,
		const Segment_2& r, RT a[], RT b[], RT c[]) const
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
	       const Segment_2& r, RT a[], RT b[], RT c[]) const 
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

    RT d[3];
    for (int i = 0; i < 3; i++) {
      d[i] = CGAL::square(l[i].a()) + CGAL::square(l[i].b());
    }

    RT z[3];
    for (int i = 0; i < 3; i++) {
      z[i] = l[(i+1)%3].a() * l[(i+2)%3].b()
	- l[(i+2)%3].a() * l[(i+1)%3].b();
    }

    
    Sqrt_1 Zero(RT(0), RT(0), d[0]);
    Sqrt_1 sqrt_D0(RT(0), RT(1), d[0]);

    Sqrt_1 D1 = d[1] + Zero;
    Sqrt_1 D2 = d[2] + Zero;

    Sqrt_3 vz(z[0] * sqrt_D0, z[1] + Zero, z[2] + Zero, Zero, D1, D2);

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

    vz = Sqrt_3(z[0] * sqrt_D0, z[1] + Zero, z[2] + Zero, Zero, D1, D2);

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

    RT x[3], y[3], w[3];
    for (int i = 0; i < 3; i++) {
      x[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      y[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      w[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);
    }

    Sqrt_3 vx, vy, vw;

    vx = Sqrt_3(x[0] * sqrt_D0, x[1] + Zero, x[2] + Zero, Zero, D1, D2);
    vy = Sqrt_3(y[0] * sqrt_D0, y[1] + Zero, y[2] + Zero, Zero, D1, D2);
    vw = Sqrt_3(w[0] * sqrt_D0, w[1] + Zero, w[2] + Zero, Zero, D1, D2);

    Sqrt_1 a1(a[(i_no+1)%3], RT(0), d[0]);
    Sqrt_1 b1(b[(i_no+1)%3], RT(0), d[0]);
    Sqrt_1 c1(c[(i_no+1)%3], RT(0), d[0]);

    Sqrt_3 dist = a1 * vx + b1 * vy + c1 * vw;

    Sign s_vw = CGAL::sign(vw);

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

    RT a[3], b[3], c[3];
    RT cx[3], cy[3], cz[3], D[3];

    orient_lines(p, q, r, a, b, c);

    for (int i = 0; i < 3; i++) {
      cx[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      cy[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      cz[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);
      D[i] = CGAL::square(a[i]) + CGAL::square(b[i]);

      Sqrt_1 Zero(RT(0), RT(0), D[0]);
      Sqrt_1 sqrt_D0(RT(0), RT(1), D[0]);
      Sqrt_1 D1 = Zero + D[1];
      Sqrt_1 D2 = Zero + D[2];

      ux = Sqrt_3(cx[0] * sqrt_D0, cx[1] + Zero, cx[2] + Zero,
		  Zero, D1, D2);
      uy = Sqrt_3(cy[0] * sqrt_D0, cy[1] + Zero, cy[2] + Zero,
		  Zero, D1, D2);
      uz = Sqrt_3(cz[0] * sqrt_D0, cz[1] + Zero, cz[2] + Zero,
		  Zero, D1, D2);      
    }
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
  //                           the orientation test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  Orientation
  orientation(const Line_2& l, PPP_Type) const
  {
    Sign s_uz = CGAL::sign(uz_ppp);
    Sign s_l =
      CGAL::sign(l.a() * ux_ppp + l.b() * uy_ppp + l.c() * uz_ppp);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }

  
  //--------------------------------------------------------------------------

  Orientation
  orientation(const Line_2& l, PPS_Type) const
  {
    Sign s_uz = CGAL::sign(uz_pps);
    Sign s_l =
      CGAL::sign(l.a() * ux_pps + l.b() * uy_pps + l.c() * uz_pps);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
  }


  //--------------------------------------------------------------------------

  // the cases PSS and SSS are identical
  template<class Type>
  Orientation
  orientation(const Line_2& l, Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().c());

    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;
    Sqrt_1 c = l.c() + Zero;

    Sign s_uz = CGAL::sign(uz);
    Sign s_l = CGAL::sign(a * ux + b * uy + c * uz);

    Sign s = Sign(s_uz * s_l);

    if ( s == ZERO ) { return COLLINEAR; }
    return ( s == POSITIVE ) ? LEFT_TURN : RIGHT_TURN;
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

  
  Sign incircle(const Point_2& t, PPS_Type type) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    Sqrt_1 vx = ux_pps - p_ref().x() * uz_pps;
    Sqrt_1 vy = uy_pps - p_ref().y() * uz_pps;

    Sqrt_1 Rs = CGAL::square(vx) + CGAL::square(vy);

    Sqrt_1 Rs1 = CGAL::square(ux_pps - t.x() * uz_pps)
      + CGAL::square(uy_pps - t.y() * uz_pps);

    return CGAL::sign(Rs1 - Rs);
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Point_2& t, PSS_Type type) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    Sqrt_1 Zero(RT(0), RT(0), ux.a().c());

    Sqrt_1 xref = p_ref().x() + Zero;
    Sqrt_1 yref = p_ref().y() + Zero;

    Sqrt_3 vx = ux - xref * uz;
    Sqrt_3 vy = uy - yref * uz;

    Sqrt_3 Rs = CGAL::square(vx) + CGAL::square(vy);

    Sqrt_1 tx = t.x() + Zero;
    Sqrt_1 ty = t.y() + Zero;

    Sqrt_3 Rs1 =
      CGAL::square(ux - tx * uz) + CGAL::square(uy - ty * uz);

    //    Sign s_vz = CGAL::sign(vz);
    Sign s_Q = CGAL::sign(Rs1 - Rs);

    //    return Sign(s_vz * s_Q);
    return s_Q;
  }


  //--------------------------------------------------------------------------

  Sign incircle(const Point_2& t, SSS_Type type) const
  {
    bool use_result(false);
    Sign s = check_easy_degeneracies(t, type, use_result);
    if ( use_result ) { return s; }

    Sqrt_1 Zero(RT(0), RT(0), ux.a().c());

    RT a1, b1, c1;
    compute_supporting_line(p_.segment(), a1, b1, c1);

    RT ns = CGAL::square(a1) + CGAL::square(b1);

    Sqrt_1 a = a1 + Zero;
    Sqrt_1 b = b1 + Zero;
    Sqrt_1 c = c1 + Zero;

    Sqrt_1 Ns = ns + Zero;
    Sqrt_3 Ls = CGAL::square(a * ux + b * uy + c * uz);

    Sqrt_1 tx = t.x() + Zero;
    Sqrt_1 ty = t.y() + Zero;

    Sqrt_3 R1s = CGAL::square(ux - tx * uz)
      + CGAL::square(uy - ty * uz);

    return CGAL::sign(R1s * Ns - Ls);
  }



  //--------------------------------------------------------------------------
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
  oriented_side(const Line_2& l, const Point_2& p, PPP_Type) const
  {
    Sign s_uz = CGAL::sign(uz_ppp);

    RT px = uz_ppp * p.x() - ux_ppp;
    RT py = uz_ppp * p.y() - uy_ppp;

    Sign s1 = CGAL::sign(l.b() * px - l.a() * py);

    Sign s = Sign(s_uz * s1);

    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }

  Oriented_side
  oriented_side(const Line_2& l, const Point_2& p, PPS_Type) const
  {
    Sqrt_1 dx = ux_pps - uz_pps * p.x();
    Sqrt_1 dy = uy_pps - uz_pps * p.y();

    Sign s = Sign(CGAL::sign(uz_pps) *
		  CGAL::sign(dy * l.a() - dx * l.b())
		  );

    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }

  // the cases PSS and SSS are identical
  template<class Type>
  Oriented_side
  oriented_side(const Line_2& l, const Point_2& p, Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().c());
    Sqrt_1 px = p.x() + Zero;
    Sqrt_1 py = p.y() + Zero;

    Sqrt_3 dx = ux - px * uz;
    Sqrt_3 dy = uy - py * uz;

    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;

    Sign s = Sign(CGAL::sign(uz) *
		  CGAL::sign(a * dy - b * dx)
		  );


    if ( s == POSITIVE ) { return ON_POSITIVE_SIDE; }
    if ( s == NEGATIVE ) { return ON_NEGATIVE_SIDE; }
    return ON_ORIENTED_BOUNDARY;
  }

  //--------------------------------------------------------------------------

  
  Sign incircle(const Line_2& l, PPP_Type) const
  {
    RT a1 = CGAL::square(l.a()) + CGAL::square(l.b());
    RT a2 = CGAL::square(ux_ppp - p_ref().x() * uz_ppp) +
      CGAL::square(uy_ppp - p_ref().y() * uz_ppp);

    RT a3 =
      CGAL::square(l.a() * ux_ppp + l.b() * uy_ppp + l.c() * uz_ppp);

    Comparison_result cr = CGAL::compare(a3, a1 * a2);

    if ( cr == LARGER ) { return POSITIVE; }
    if ( cr == SMALLER ) { return NEGATIVE; }
    return ZERO;
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, PPS_Type) const
  {
    Sqrt_1 vx = ux_pps - p_ref().x() * uz_pps;
    Sqrt_1 vy = uy_pps - p_ref().y() * uz_pps;

    Sqrt_1 Rs = CGAL::square(vx) + CGAL::square(vy);
    
    RT Ns = CGAL::square(l.a()) + CGAL::square(l.b());

    Sqrt_1 Ls =
      CGAL::square(l.a() * ux_pps + l.b() * uy_pps + l.c() * uz_pps);

    return CGAL::sign(Ls - Rs * Ns);
  }


  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, PSS_Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().c());

    Sqrt_3 vx = ux - (p_ref().x() + Zero) * uz;
    Sqrt_3 vy = uy - (p_ref().y() + Zero) * uz;

    Sqrt_3 Rs = CGAL::square(vx) + CGAL::square(vy);


    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;
    Sqrt_1 c = l.c() + Zero;

    Sqrt_1 Ns = CGAL::square(a) + CGAL::square(b);

    Sqrt_3 Ls = CGAL::square(a * ux + b * uy + c * uz);

    return CGAL::sign(Ls - Rs * Ns);
  }


  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, SSS_Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().c());

    RT a1, b1, c1;
    compute_supporting_line(p_.segment(), a1, b1, c1);

    Sqrt_1 a = a1 + Zero;
    Sqrt_1 b = b1 + Zero;
    Sqrt_1 c = c1 + Zero;

    Sqrt_1 Ns = CGAL::square(a) + CGAL::square(b);

    Sqrt_1 la = l.a() + Zero;
    Sqrt_1 lb = l.b() + Zero;
    Sqrt_1 lc = l.c() + Zero;

    Sqrt_1 Ns1 = CGAL::square(la) + CGAL::square(lb);

    Sqrt_3 Ls = CGAL::square(a * ux + b * uy + c * uz);

    Sqrt_3 Ls1 =
      CGAL::square(la * ux + lb * uy + lc * uz);

    return CGAL::sign(Ls1 * Ns - Ls * Ns1);
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------



  template<class Type>
  Sign incircle(const Segment_2& t, Type type) const
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
    Sign sl = incircle(l, type);

    // this is the old code
    if ( sl == POSITIVE ) { return sl; }

    Oriented_side os1 = oriented_side(l, t.source(), type);
    Oriented_side os2 = oriented_side(l, t.target(), type);

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
  bool is_same_point(const Point_2& p, Method_tag) const
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
    if ( v_type == PPP ) { return ux_ppp; }
    if ( v_type == PPS ) { return to_ft(ux_pps); }
    return to_ft(ux);
  }

  FT hy() const {
    if ( v_type == PPP ) { return uy_ppp; }
    if ( v_type == PPS ) { return to_ft(uy_pps); }
    return to_ft(uy);
  }

  FT hw() const {
    if ( v_type == PPP ) { return uz_ppp; }
    if ( v_type == PPS ) { return to_ft(uz_pps); }
    return to_ft(uz);
  }

  FT squared_radius() const {
    switch (v_type) {
    case PPP:    case PPS:    case PSS:
      {
	Point_2 p_ref;
	if ( p_.is_point() ) {
	  p_ref = p_.point();
	} else if ( q_.is_point() ) {
	  p_ref = q_.point();
	} else {
	  CGAL_assertion( r_.is_point() );
	  p_ref = r_.point();
	}

	FT dx2 = CGAL::square(x() - p_ref.x());
	FT dy2 = CGAL::square(y() - p_ref.y());
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
  Svd_voronoi_vertex_ring_C2(const Site_2& p,
			     const Site_2& q,
			     const Site_2& r)
    : p_(p), q_(q), r_(r)
  {
    compute_vertex(p, q, r);
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Site_2& t) const 
  {
    Sign s;

    if ( t.is_point() ) {
      s = incircle(t.point());
    } else {
      s = incircle(t.segment());
    }

    return s;
  }

  //--------------------------------------------------------------------------


  Orientation orientation(const Line_2& l) const 
  {
    Orientation o(COLLINEAR);
    switch ( v_type ) {
    case PPP:
      o = orientation(l, PPP_Type());
      break;
    case PPS:
      o = orientation(l, PPS_Type());
      break;
    case PSS:
      o = orientation(l, PSS_Type());
      break;
    case SSS:
      o = orientation(l, SSS_Type());
      break;
    }

    return o;
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

  // the case ppp
  RT ux_ppp, uy_ppp, uz_ppp;

  // the case pps
  Sqrt_1 ux_pps, uy_pps, uz_pps;

  // the case pss and sss
  Sqrt_3 ux, uy, uz;
};








CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_RING_C2_H

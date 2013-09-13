// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_RING_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_RING_C2_H


#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Kernel/global_functions_2.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {


template<class K>
class Voronoi_vertex_ring_C2
  : public Basic_predicates_C2<K>
{
public:
  typedef Basic_predicates_C2<K> Base;

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

  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Oriented_side       Oriented_side;
  typedef typename Base::Sign                Sign;

  using Base::compute_supporting_line;
  using Base::oriented_side_of_line;
  using Base::opposite_line;
  using Base::to_ft;

private:
  typedef Are_same_points_C2<K>     Are_same_points_2;
  typedef Are_same_segments_C2<K>   Are_same_segments_2;

  typedef typename K::Intersections_tag ITag;

  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

private:
  //--------------------------------------------------------------------------

  void
  compute_ppp(const Site_2& sp, const Site_2& sq, const Site_2& sr)
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
		       sr.is_point() );

    Point_2 p = sp.point(), q = sq.point(), r = sr.point();

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
  compute_pss(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_point() && q.is_segment() &&
		       r.is_segment() );

    v_type = PSS;

    bool pq =
      same_points(p, q.source_site()) || same_points(p, q.target_site());
    bool pr =
      same_points(p, r.source_site()) || same_points(p, r.target_site());

    Point_2 pp = p.point();

    if ( pq && pr ) {
      Sqrt_1 One(RT(1), RT(0), RT(0));

      ux = Sqrt_3(pp.x() * One);
      uy = Sqrt_3(pp.y() * One);
      uz = Sqrt_3(One);
      return;
    }

    

    RT a1, b1, c1, a2, b2, c2;
    compute_supporting_line(q.supporting_site(), a1, b1, c1);
    compute_supporting_line(r.supporting_site(), a2, b2, c2);

    RT c1_ = a1 * pp.x() + b1 * pp.y() + c1;
    RT c2_ = a2 * pp.x() + b2 * pp.y() + c2;

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
      if ( same_points(p, q.target_site()) ) {
	a1 = -a1;  b1 = -b1;  c1_ = -c1_;
      }
    }

    if ( sgn_c2_ == NEGATIVE ) {
      a2 = -a2;  b2 = -b2;  c2_ = -c2_;
    } else if ( sgn_c2_ == ZERO ) {
      CGAL_assertion( pr );
      if ( same_points(p, r.source_site()) ) {
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

      ux = Sqrt_3(J + pp.x() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uy = Sqrt_3(I + pp.y() * vz,
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

      ux = Sqrt_3(J + pp.x() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uy = Sqrt_3(I + pp.y() * vz,
		  Zero, Zero, Zero, Zero, Zero);
      uz = Sqrt_3(vz, Zero, Zero, Zero, Zero, Zero);
    } else {
      Line_2 lq(a1, b1, c1_);
      Line_2 lr(a2, b2, c2_);
      compute_pll(pp, lq, lr);
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
      s_sigma = s1 * CGAL::sign(e);
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
      s_rho = s1 * CGAL::sign(e);
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
  compute_pps(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_point() && q.is_point() &&
		       r.is_segment() );

    v_type = PPS;

    RT a, b, c;
    compute_supporting_line(r.supporting_site(), a, b, c);

    Point_2 pp = p.point(), qq = q.point();

    RT c_ = a * pp.x() + b * pp.y() + c;
    RT cq_ = a * qq.x() + b * qq.y() + c;


    if ( same_points(p, r.source_site()) ||
	 same_points(p, r.target_site()) ) {
      c_ = RT(0);
    }
    if ( same_points(q, r.source_site()) ||
	 same_points(q, r.target_site()) ) {
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

    RT x_ = qq.x() - pp.x();
    RT y_ = qq.y() - pp.y();
    RT n_ = CGAL::square(x_) + CGAL::square(y_);


    Comparison_result res = CGAL::compare( c_, cq_ );

    if ( res == EQUAL ) {
      RT e1 = CGAL::square(c_);
      RT J = nl * (a * n_ + RT(4) * c_ * x_) - RT(4) * a * e1;
      RT I = nl * (b * n_ + RT(4) * c_ * y_) - RT(4) * b * e1;
      RT X = RT(8) * nl * c_;

      ux_pps = Sqrt_1(J + pp.x() * X);
      uy_pps = Sqrt_1(I + pp.y() * X);
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

    ux_pps = Sqrt_1(J + pp.x() * X, RT(-2) * y_, S);
    uy_pps = Sqrt_1(I + pp.y() * X, RT( 2) * x_, S);
    uz_pps = Sqrt_1(X, RT(0), S);
  }

  //--------------------------------------------------------------------------

  bool check_if_exact(const Site_2& , unsigned int ,
		      const Tag_false&) const
  {
    return true;
  }

  bool check_if_exact(const Site_2& s, unsigned int i,
		      const Tag_true&) const
  {
    return s.is_input(i);
  }

  // determines of the segment s is on the positive halfspace as
  // defined by the supporting line of the segment supp; the line l
  // is supposed to be the supporting line of the segment supp and we
  // pass it so that we do not have to recompute it
  bool
  is_on_positive_halfspace(const Site_2& supp,
			   const Site_2& s, const Line_2& l) const
  {
    CGAL_precondition( supp.is_segment() && s.is_segment() );
    
    if ( same_segments(supp.supporting_site(),
		       s.supporting_site()) ) {
      return false;
    }

    if ( same_points(supp.source_site(), s.source_site()) ||
	 same_points(supp.target_site(), s.source_site()) ) {
      return oriented_side_of_line(l, s.target()) == ON_POSITIVE_SIDE;
    }

    if ( same_points(supp.source_site(), s.target_site()) ||
	 same_points(supp.target_site(), s.target_site()) ) {
      return oriented_side_of_line(l, s.source()) == ON_POSITIVE_SIDE;
    }

    ITag itag;

    if ( !check_if_exact(s, 0, itag) &&
	 same_segments(supp.supporting_site(),
		       s.crossing_site(0)) ) {
      return oriented_side_of_line(l, s.target()) == ON_POSITIVE_SIDE;
    }

    if ( !check_if_exact(s, 1, itag) &&
	 same_segments(supp.supporting_site(),
		       s.crossing_site(1)) ) {
      return oriented_side_of_line(l, s.source()) == ON_POSITIVE_SIDE;
    }

    return Base::is_on_positive_halfspace(l, s.segment());
  }

  //--------------------------------------------------------------------------

  void
  orient_lines(const Site_2& p, const Site_2& q,
	       const Site_2& r, RT a[], RT b[], RT c[]) const 
  {
    CGAL_precondition( p.is_segment() && q.is_segment() &&
		       r.is_segment() );

    Line_2 l[3];
    l[0] = compute_supporting_line(p.supporting_site());
    l[1] = compute_supporting_line(q.supporting_site());
    l[2] = compute_supporting_line(r.supporting_site());
    
    bool is_oriented[3] = {false, false, false};

    if ( is_on_positive_halfspace(p, q, l[0]) ||
	 is_on_positive_halfspace(p, r, l[0]) ) {
      is_oriented[0] = true;
    } else {
      l[0] = opposite_line(l[0]);
      if ( is_on_positive_halfspace(p, q, l[0]) ||
	   is_on_positive_halfspace(p, r, l[0]) ) {
	is_oriented[0] = true;
      } else {
	l[0] = opposite_line(l[0]);
      }
    }

    if ( is_on_positive_halfspace(q, p, l[1]) ||
	 is_on_positive_halfspace(q, r, l[1]) ) {
      is_oriented[1] = true;
    } else {
       l[1] = opposite_line(l[1]);
      if ( is_on_positive_halfspace(q, p, l[1]) ||
	   is_on_positive_halfspace(q, r, l[1]) ) {
	is_oriented[1] = true;
      } else {
	l[1] = opposite_line(l[1]);
      }
    }

    if ( is_on_positive_halfspace(r, p, l[2]) ||
	 is_on_positive_halfspace(r, q, l[2]) ) {
      is_oriented[2] = true;
    } else {
      l[2] = opposite_line(l[2]);
      if ( is_on_positive_halfspace(r, p, l[2]) ||
	   is_on_positive_halfspace(r, q, l[2]) ) {
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

    Sign sgn_dist = s_vw * CGAL::sign(dist);

    CGAL_assertion( sgn_dist != ZERO );

    if ( sgn_dist == NEGATIVE ) {
      a[i_no] = -a[i_no];
      b[i_no] = -b[i_no];
      c[i_no] = -c[i_no];
    }
  }


  void
  compute_sss(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_segment() && q.is_segment() &&
		       r.is_segment() );

    v_type = SSS;

    RT a[3], b[3], c[3];
    RT cx[3], cy[3], cz[3], D[3];

    orient_lines(p, q, r, a, b, c);

    for (int i = 0; i < 3; i++) {
      cx[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      cy[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      cz[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);
      D[i] = CGAL::square(a[i]) + CGAL::square(b[i]);
    }

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

  //--------------------------------------------------------------------------

  void
  compute_vertex(const Site_2& s1, const Site_2& s2, const Site_2& s3)
  {
    if ( s1.is_point() && s2.is_point() && s3.is_point() ) {
      compute_ppp(s1, s2, s3);

    } else if ( s1.is_segment() && s2.is_point() && s3.is_point() ) {
      compute_vertex(s2, s3, s1);
      pps_idx = 1;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
      pps_idx = 2;

    } else if ( s1.is_point() && s2.is_point() && s3.is_segment() ) {
      compute_pps(s1, s2, s3);
      pps_idx = 0;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_segment() ) {
      compute_pss(s1, s2, s3);
    } else if ( s1.is_segment() && s2.is_point() && s3.is_segment() ) {
      compute_vertex(s2, s3, s1);
    } else if ( s1.is_segment() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
    } else {
      compute_sss(s1, s2, s3);
    }

  }

  //--------------------------------------------------------------------------

  bool is_endpoint_of(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    return ( same_points(p, s.source_site()) ||
	     same_points(p, s.target_site()) );
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

    return s_uz * s_l;
  }

  
  //--------------------------------------------------------------------------

  Orientation
  orientation(const Line_2& l, PPS_Type) const
  {
    Sign s_uz = CGAL::sign(uz_pps);
    Sign s_l = CGAL::sign(l.a() * ux_pps + l.b() * uy_pps + l.c() * uz_pps);

    return s_uz * s_l;
  }


  //--------------------------------------------------------------------------

  // the cases PSS and SSS are identical
  template<class Type>
  Orientation
  orientation(const Line_2& l, Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().root());

    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;
    Sqrt_1 c = l.c() + Zero;

    Sign s_uz = CGAL::sign(uz);
    Sign s_l = CGAL::sign(a * ux + b * uy + c * uz);

    return s_uz * s_l;
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

  Sign check_easy_degeneracies(const Site_2& t, PPS_Type,
			       bool& use_result) const
  {
    CGAL_precondition( t.is_point() );

    use_result = false;
    if (  ( p_.is_point() && same_points(p_, t) ) ||
	  ( q_.is_point() && same_points(q_, t) ) ||
	  ( r_.is_point() && same_points(r_, t) )  ) {
      use_result = true;
      return ZERO;
    }

    if (  ( p_.is_segment() && is_endpoint_of(t, p_) ) || 
	  ( q_.is_segment() && is_endpoint_of(t, q_) ) ||
	  ( r_.is_segment() && is_endpoint_of(t, r_) )  ) {
      use_result = true;
      return POSITIVE;
    }

    return ZERO;
  }

  inline
  Sign check_easy_degeneracies(const Site_2& t, PSS_Type,
			       bool& use_result) const
  {
    CGAL_precondition( t.is_point() );
    return check_easy_degeneracies(t, PPS_Type(), use_result);
  }

  inline
  Sign check_easy_degeneracies(const Site_2& t, SSS_Type,
			       bool& use_result) const
  {
    CGAL_USE(t);
    CGAL_precondition( t.is_point() );
    use_result = false;
    // ADD THE CASES WHERE t IS AN ENDPOINT OF ONE OF THE SEGMENTS
    return ZERO;
  }

  //--------------------------------------------------------------------------

  template<class Type>
  inline  
  Sign incircle_p(const Site_2& st, Type type) const
  {
    CGAL_precondition( st.is_point() );

    bool use_result(false);
    Sign s = check_easy_degeneracies(st, type, use_result);
    if ( use_result ) { return s; }

    return incircle_p_no_easy(st, type);
  }

  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& st, PPP_Type) const
  {
    CGAL_precondition( st.is_point() );

    Point_2 t = st.point();

    return - side_of_oriented_circle(p_.point(), q_.point(), r_.point(), t);
  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, PPS_Type ) const
  {
    CGAL_precondition( st.is_point() );

    Point_2 t = st.point();

    Point_2 pref = p_ref().point();

    Sqrt_1 vx = ux_pps - pref.x() * uz_pps;
    Sqrt_1 vy = uy_pps - pref.y() * uz_pps;

    Sqrt_1 Rs = CGAL::square(vx) + CGAL::square(vy);

    Sqrt_1 Rs1 = CGAL::square(ux_pps - t.x() * uz_pps)
      + CGAL::square(uy_pps - t.y() * uz_pps);

    return CGAL::sign(Rs1 - Rs);
  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, PSS_Type ) const
  {
    CGAL_precondition( st.is_point() );
    Point_2 t = st.point();

    Sqrt_1 Zero(RT(0), RT(0), ux.a().root());

    Point_2 pref = p_ref().point();

    Sqrt_1 xref = pref.x() + Zero;
    Sqrt_1 yref = pref.y() + Zero;

    Sqrt_3 vx = ux - xref * uz;
    Sqrt_3 vy = uy - yref * uz;

    Sqrt_3 Rs = CGAL::square(vx) + CGAL::square(vy);

    Sqrt_1 tx = t.x() + Zero;
    Sqrt_1 ty = t.y() + Zero;

    Sqrt_3 Rs1 =
      CGAL::square(ux - tx * uz) + CGAL::square(uy - ty * uz);

    Sign s_Q = CGAL::sign(Rs1 - Rs);

    return s_Q;
  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, SSS_Type ) const
  {
    CGAL_precondition( st.is_point() );

    Point_2 t = st.point();

    Sqrt_1 Zero(RT(0), RT(0), ux.a().root());

    RT a1, b1, c1;
    compute_supporting_line(p_.supporting_site(), a1, b1, c1);

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

  Sign incircle_p(const Site_2& t) const 
  {
    if ( is_degenerate_Voronoi_circle() ) {
      return POSITIVE;
    }

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_p(t, PPP_Type());
      break;
    case PPS:
      s = incircle_p(t, PPS_Type());
      break;
    case PSS:
      s = incircle_p(t, PSS_Type());
      break;
    case SSS:
      s = incircle_p(t, SSS_Type());
      break;
    }

    return s;
  }

  Sign incircle_p_no_easy(const Site_2& t) const 
  {
    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_p(t, PPP_Type());
      break;
    case PPS:
      s = incircle_p_no_easy(t, PPS_Type());
      break;
    case PSS:
      s = incircle_p_no_easy(t, PSS_Type());
      break;
    case SSS:
      s = incircle_p_no_easy(t, SSS_Type());
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

    return s_uz * s1;
  }

  Oriented_side
  oriented_side(const Line_2& l, const Point_2& p, PPS_Type) const
  {
    Sqrt_1 dx = ux_pps - uz_pps * p.x();
    Sqrt_1 dy = uy_pps - uz_pps * p.y();

    return CGAL::sign(uz_pps) * CGAL::sign(dy * l.a() - dx * l.b());
  }

  // the cases PSS and SSS are identical
  template<class Type>
  Oriented_side
  oriented_side(const Line_2& l, const Point_2& p, Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().root());
    Sqrt_1 px = p.x() + Zero;
    Sqrt_1 py = p.y() + Zero;

    Sqrt_3 dx = ux - px * uz;
    Sqrt_3 dy = uy - py * uz;

    Sqrt_1 a = l.a() + Zero;
    Sqrt_1 b = l.b() + Zero;

    return CGAL::sign(uz) * CGAL::sign(a * dy - b * dx);
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  
  Sign incircle(const Line_2& l, PPP_Type) const
  {
    Point_2 pref = p_ref().point();

    RT a1 = CGAL::square(l.a()) + CGAL::square(l.b());
    RT a2 = CGAL::square(ux_ppp - pref.x() * uz_ppp) +
      CGAL::square(uy_ppp - pref.y() * uz_ppp);

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
    Point_2 pref = p_ref().point();

    Sqrt_1 vx = ux_pps - pref.x() * uz_pps;
    Sqrt_1 vy = uy_pps - pref.y() * uz_pps;

    Sqrt_1 Rs = CGAL::square(vx) + CGAL::square(vy);
    
    RT Ns = CGAL::square(l.a()) + CGAL::square(l.b());

    Sqrt_1 Ls =
      CGAL::square(l.a() * ux_pps + l.b() * uy_pps + l.c() * uz_pps);

    return CGAL::sign(Ls - Rs * Ns);
  }


  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, PSS_Type) const
  {
    Sqrt_1 Zero(RT(0), RT(0), ux.a().root());

    Point_2 pref = p_ref().point();

    Sqrt_3 vx = ux - (pref.x() + Zero) * uz;
    Sqrt_3 vy = uy - (pref.y() + Zero) * uz;

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
    Sqrt_1 Zero(RT(0), RT(0), ux.a().root());

    RT a1, b1, c1;
    compute_supporting_line(p_.supporting_site(), a1, b1, c1);

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
  Sign incircle_s(const Site_2& t, Type type) const
  {
    CGAL_precondition( t.is_segment() );

    if ( v_type == PPP || v_type == PPS ) {
      if (  p_.is_point() && q_.is_point() &&
	    is_endpoint_of(p_, t) && is_endpoint_of(q_, t)  ) {
	return NEGATIVE;
      }

      if (  p_.is_point() && r_.is_point() &&
	    is_endpoint_of(p_, t) && is_endpoint_of(r_, t)  ){
	return NEGATIVE;
      }

      if (  q_.is_point() && r_.is_point() &&
	    is_endpoint_of(q_, t) && is_endpoint_of(r_, t)  ){
	return NEGATIVE;
      }
    }

    if ( v_type == PSS ) {
      if ( p_.is_segment() &&
	   same_segments(p_.supporting_site(),
			 t.supporting_site()) ) {
	return POSITIVE;
      }
      if ( q_.is_segment() &&
	   same_segments(q_.supporting_site(),
			 t.supporting_site()) ) {
	return POSITIVE;
      }
      if ( r_.is_segment() &&
	   same_segments(r_.supporting_site(),
			 t.supporting_site()) ) {
	return POSITIVE;
      }
    }

    return incircle_s_no_easy(t, type);
  }

  template<class Type>
  Sign incircle_s_no_easy(const Site_2& t, Type type) const
  {
    Sign d1, d2;
    if (  ( p_.is_point() && same_points(p_, t.source_site()) ) ||
	  ( q_.is_point() && same_points(q_, t.source_site()) ) ||
	  ( r_.is_point() && same_points(r_, t.source_site()) )  ) {
      d1 = ZERO;
    } else {
      d1 = incircle_p(t.source_site());
    }
    if ( d1 == NEGATIVE ) { return NEGATIVE; }

    if (  ( p_.is_point() && same_points(p_, t.target_site()) ) ||
	  ( q_.is_point() && same_points(q_, t.target_site()) ) ||
	  ( r_.is_point() && same_points(r_, t.target_site()) )  ) {
      d2 = ZERO;
    } else {
      d2 = incircle_p(t.target_site());
    }
    if ( d2 == NEGATIVE ) { return NEGATIVE; }

    Line_2 l = compute_supporting_line(t.supporting_site());
    Sign sl = incircle(l, type);

    if ( sl == POSITIVE ) { return sl; }

    if ( sl == ZERO && (d1 == ZERO || d2 == ZERO) ) { return ZERO; }

    Oriented_side os1 = oriented_side(l, t.source(), type);
    Oriented_side os2 = oriented_side(l, t.target(), type);

    if ( sl == ZERO ) {
      if (os1 == ON_ORIENTED_BOUNDARY || os2 == ON_ORIENTED_BOUNDARY) {
	return ZERO;
      }
      return ( os1 == os2 ) ? POSITIVE : ZERO;
    }

    return (os1 == os2) ? POSITIVE : NEGATIVE;
  }

  //--------------------------------------------------------------------------

  

  Sign incircle_s(const Site_2& t) const 
  {
    CGAL_precondition( t.is_segment() );

    if ( is_degenerate_Voronoi_circle() ) {
      // case 1: the new segment is not adjacent to the center of the
      //         degenerate Voronoi circle
      if (  !same_points( p_ref(), t.source_site() ) &&
	    !same_points( p_ref(), t.target_site() )  ) {
	return POSITIVE;
      }

      CGAL_assertion( v_type == PSS );

      if ( p_.is_segment() &&
	   same_segments(p_.supporting_site(),
			 t.supporting_site()) ) {
	return ZERO;
      }

      if ( q_.is_segment() &&
	   same_segments(q_.supporting_site(),
			 t.supporting_site()) ) {
	return ZERO;
      }

      if ( r_.is_segment() &&
	   same_segments(r_.supporting_site(),
			 t.supporting_site()) ) {
	return ZERO;
      }

      Site_2 pr;
      Site_2 sp, sq;
      if ( p_.is_point() ) {
	CGAL_assertion( q_.is_segment() && r_.is_segment() );
	pr = p_;
	sp = q_;
	sq = r_;
      } else if ( q_.is_point() ) {
	CGAL_assertion( r_.is_segment() && p_.is_segment() );
	pr = q_;
	sp = r_;
	sq = p_;
      } else {
	CGAL_assertion( p_.is_segment() && q_.is_segment() );
	pr = r_;
	sp = p_;
	sq = q_;
      }

      Point_2 pq = sq.source(), pp = sp.source(), pt = t.source();

      if ( same_points(sp.source_site(), pr) ) { pp = sp.target(); }
      if ( same_points(sq.source_site(), pr) ) { pq = sq.target(); }
      if ( same_points( t.source_site(), pr) ) { pt =  t.target(); }

      Point_2 pr_ = pr.point();

      if ( CGAL::orientation(pr_, pp, pt) == LEFT_TURN &&
	   CGAL::orientation(pr_, pq, pt) == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return ZERO;
    } // if ( is_degenerate_Voronoi_circle() )

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_s(t, PPP_Type());
      break;
    case PPS:
      s = incircle_s(t, PPS_Type());
      break;
    case PSS:
      s = incircle_s(t, PSS_Type());
      break;
    case SSS:
      s = incircle_s(t, SSS_Type());
      break;
    }

    return s;
  }

  Sign incircle_s_no_easy(const Site_2& t) const
  {
    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_s_no_easy(t, PPP_Type());
      break;
    case PPS:
      s = incircle_s_no_easy(t, PPS_Type());
      break;
    case PSS:
      s = incircle_s_no_easy(t, PSS_Type());
      break;
    case SSS:
      s = incircle_s_no_easy(t, SSS_Type());
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
      return ( is_endpoint_of(p_, q_) && is_endpoint_of(p_, r_) );
    } else if ( q_.is_point() ) {
      return ( is_endpoint_of(q_, p_) && is_endpoint_of(q_, r_) );
    } else {
      CGAL_assertion( r_.is_point() );
      return ( is_endpoint_of(r_, p_) && is_endpoint_of(r_, q_) );
    }
  }


  //--------------------------------------------------------------------------

private:

  //--------------------------------------------------------------------------
  //  the reference point (valid if v_type != SSS)
  //--------------------------------------------------------------------------

  Site_2 p_ref() const
  {
    CGAL_precondition ( v_type != SSS );

    if ( v_type == PPS ) {
      if ( pps_idx == 0 ) { return p_; }
      if ( pps_idx == 1 ) { return q_; }
      return r_;
    }

    if ( p_.is_point() ) {
      return p_;
    } else if ( q_.is_point() ) {
      return q_;
    } else {
      CGAL_assertion( r_.is_point() );
      return r_;
    }
  }


public:
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                           access methods
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  inline FT x(Integral_domain_without_division_tag) const {
    return CGAL::to_double(hx()) / CGAL::to_double(hw());
  }
  inline FT y(Integral_domain_without_division_tag) const {
    return CGAL::to_double(hy()) / CGAL::to_double(hw());
  }

  inline FT x(Field_tag) const { return hx() / hw(); }
  inline FT y(Field_tag) const { return hy() / hw(); }

  inline FT x() const {
      typedef Algebraic_structure_traits<FT> AST;
      return x(typename AST::Algebraic_category());
  }
    
  inline FT y() const {
      typedef Algebraic_structure_traits<FT> AST;
      return y(typename AST::Algebraic_category());
  }

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
	Point_2 pref = p_ref().point();
	FT dx2 = CGAL::square(x() - pref.x());
	FT dy2 = CGAL::square(y() - pref.y());
	return dx2 + dy2;
      }
      break;
    case SSS:
      {
	Line_2 l = compute_supporting_line(p_.supporting_site());
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
      return degenerate_point();
    }
    
    return Point_2(x(), y());
  }


  Point_2 degenerate_point() const
  {
    CGAL_precondition( is_degenerate_Voronoi_circle() );
    return p_ref().point();
  }


  typename K::Circle_2 circle() const
  {
    typedef typename K::Circle_2  K_Circle_2;
    return K_Circle_2(point(), squared_radius());
  }

  vertex_t type() const { return v_type; }

public:
  Voronoi_vertex_ring_C2(const Site_2& p,
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
      s = incircle_p(t);
    } else {
      s = incircle_s(t);
    }

    return s;
  }

  Sign incircle_no_easy(const Site_2& t) const 
  {
    Sign s;

    if ( t.is_point() ) {
      s = incircle_p_no_easy(t);
    } else {
      s = incircle_s_no_easy(t);
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


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_RING_C2_H

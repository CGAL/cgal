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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_SQRT_FIELD_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_SQRT_FIELD_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>




#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {


template<class K>
class Voronoi_vertex_sqrt_field_C2
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

  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Oriented_side       Oriented_side;
  typedef typename Base::Sign                Sign;
  typedef typename Base::Compute_scalar_product_2 Compute_scalar_product_2;

private:
  typedef Are_same_points_C2<K>    Are_same_points_2;
  typedef Are_same_segments_C2<K>  Are_same_segments_2;

  typedef typename K::Intersections_tag ITag;

  Are_same_points_2     same_points;
  Are_same_segments_2   same_segments;


private:
  //--------------------------------------------------------------------------

  void
  compute_ppp(const Site_2& sp, const Site_2& sq, const Site_2& sr)
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
		       sr.is_point() );

    Point_2 p = sp.point(), q = sq.point(), r = sr.point();

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
  compute_pss(const Site_2& p, const Site_2& q, const Site_2& r)
  {
#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      std::ofstream ofs("vv-failure-log.cin", std::ios_base::app);
      ofs.precision(16);
      ofs << p << std::endl;
      ofs << q << std::endl;
      ofs << r << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

    CGAL_precondition( p.is_point() && q.is_segment() &&
		       r.is_segment() );

    bool pq =
      same_points(p, q.source_site()) || same_points(p, q.target_site());
    bool pr =
      same_points(p, r.source_site()) || same_points(p, r.target_site());

    Point_2 pp = p.point();

    if ( pq && pr ) {
      ux = pp.x();
      uy = pp.y();
      uz = FT(1);
      return;
    }


    FT a1, b1, c1, a2, b2, c2;
    compute_supporting_line(q.supporting_site(), a1, b1, c1);
    compute_supporting_line(r.supporting_site(), a2, b2, c2);

    FT c1_ = a1 * pp.x() + b1 * pp.y() + c1;
    FT c2_ = a2 * pp.x() + b2 * pp.y() + c2;

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
      FT J = a1 * c2_;
      FT I = b1 * c2_;

      FT n1 = CGAL::square(a1) + CGAL::square(b1);
      FT n2 = CGAL::square(a2) + CGAL::square(b2);

      FT D1D2 = n1 * n2;

      uz = -a1 * a2 - b1 * b2 + CGAL::sqrt(D1D2);

      ux = J + pp.x() * uz;
      uy = I + pp.y() * uz;

    } else if ( pr ) {
      FT J = a2 * c1_;
      FT I = b2 * c1_;

      FT n1 = CGAL::square(a1) + CGAL::square(b1);
      FT n2 = CGAL::square(a2) + CGAL::square(b2);

      FT D1D2 = n1 * n2;

      uz = -a1 * a2 - b1 * b2 + CGAL::sqrt(D1D2);

      ux = J + pp.x() * uz;
      uy = I + pp.y() * uz;

    } else {
      Line_2 lq(a1, b1, c1_);
      Line_2 lr(a2, b2, c2_);
      compute_pll(pp, lq, lr);
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
  compute_pps(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_point() && q.is_point() &&
		       r.is_segment() );
    FT a, b, c;
    compute_supporting_line(r.supporting_site(), a, b, c);

    Point_2 pp = p.point(), qq = q.point();

    FT c_ = a * pp.x() + b * pp.y() + c;
    FT cq_ = a * qq.x() + b * qq.y() + c;

    if ( same_points(p, r.source_site()) ||
	 same_points(p, r.target_site()) ) {
      c_ = FT(0);
    }
    if ( same_points(q, r.source_site()) ||
	 same_points(q, r.target_site()) ) {
      cq_ = FT(0);
    }

#if 0
    // MK:: 1/3/2010
    // eliminatine the following code seems to have no effect. the
    // analysis supports this; this code is to be removed after some
    // testing.
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
#endif

    FT nl = CGAL::square(a) + CGAL::square(b);

    FT x_ = qq.x() - pp.x();
    FT y_ = qq.y() - pp.y();
    FT n_ = CGAL::square(x_) + CGAL::square(y_);

    // the following lines of code check whether the line of the two
    // points is parallel to the supporting line of the segment, and
    // parallel to the x or y axis.

    Point_2 r_src = r.source_site().point();
    Point_2 r_trg = r.target_site().point();

    bool pq_xaligned = pp.y() == qq.y();
    bool r_xaligned = r_src.y() == r_trg.y();

    bool pq_yaligned = pp.x() == qq.x();
    bool r_yaligned = r_src.x() == r_trg.x();

    bool parallel = (pq_xaligned && r_xaligned) || (pq_yaligned && r_yaligned);

    if ( parallel || c_ == cq_ ) {
      FT e1 = CGAL::square(c_);
      FT J = nl * (a * n_ + FT(4) * c_ * x_) - FT(4) * a * e1;
      FT I = nl * (b * n_ + FT(4) * c_ * y_) - FT(4) * b * e1;
      FT X = FT(8) * nl * c_;

      ux = J + pp.x() * X;
      uy = I + pp.y() * X;
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

    ux = J + pp.x() * X - FT(2) * y_ * sqrt_S;
    uy = I + pp.y() * X + FT(2) * x_ * sqrt_S;
    uz = X;
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
  orient_lines(const Site_2& sp, const Site_2& sq,
	       const Site_2& sr, FT a[], FT b[], FT c[]) const 
  {
    CGAL_precondition( sp.is_segment() && sq.is_segment() &&
		       sr.is_segment() );

    Line_2 l[3];
    l[0] = compute_supporting_line(sp.supporting_site());
    l[1] = compute_supporting_line(sq.supporting_site());
    l[2] = compute_supporting_line(sr.supporting_site());
    
    bool is_oriented[3] = {false, false, false};

    if ( is_on_positive_halfspace(sp, sq, l[0]) ||
    	 is_on_positive_halfspace(sp, sr, l[0]) ) {
      is_oriented[0] = true;
    } else {
      
      l[0] = opposite_line(l[0]);
      if ( is_on_positive_halfspace(sp, sq, l[0]) ||
      	   is_on_positive_halfspace(sp, sr, l[0]) ) {
	is_oriented[0] = true;
      } else {
	l[0] = opposite_line(l[0]);
      }
    }

    if ( is_on_positive_halfspace(sq, sp, l[1]) ||
	 is_on_positive_halfspace(sq, sr, l[1]) ) {
      is_oriented[1] = true;
    } else {
      l[1] = opposite_line(l[1]);
      if ( is_on_positive_halfspace(sq, sp, l[1]) ||
	   is_on_positive_halfspace(sq, sr, l[1]) ) {
	is_oriented[1] = true;
      } else {
	l[1] = opposite_line(l[1]);
      }
    }

    if ( is_on_positive_halfspace(sr, sp, l[2]) ||
	 is_on_positive_halfspace(sr, sq, l[2]) ) {
      is_oriented[2] = true;
    } else {
      l[2] = opposite_line(l[2]);
      if ( is_on_positive_halfspace(sr, sp, l[2]) ||
	   is_on_positive_halfspace(sr, sq, l[2]) ) {
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
    int npts = 0;
    if ( s1.is_point() ) ++npts;
    if ( s2.is_point() ) ++npts;
    if ( s3.is_point() ) ++npts;

    switch ( npts ) {
    case 0:
      v_type = SSS;
      break;
    case 1:
      v_type = PSS;
      break;
    case 2:
      v_type = PPS;
      break;
    default:
      v_type = PPP;
    }


    if ( v_type == PPP ) {
      compute_ppp(s1, s2, s3);
    } else if ( v_type == SSS ) {
      compute_sss(s1, s2, s3);
    } else if ( v_type == PPS ) {
      if ( s1.is_segment() ) {
	compute_pps(s2, s3, s1);
	pps_idx = 1;
      } else if ( s2.is_segment() ) {
	compute_pps(s3, s1, s2);
	pps_idx = 2;
      } else {
	compute_pps(s1, s2, s3);
	pps_idx = 0;
      }
    } else {
      if ( s1.is_point() ) {
	compute_pss(s1, s2, s3);
      } else if ( s2.is_point() ) {
	compute_pss(s2, s3, s1);
      } else {
	compute_pss(s3, s1, s2);
      }
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
  //                              the incircle test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a point
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------

  Sign
  check_easy_degeneracies(const Site_2& t, PPS_Type,
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

  Sign
  check_easy_degeneracies(const Site_2& t, PSS_Type,
			  bool& use_result) const
  {
    CGAL_precondition( t.is_point() );

    return check_easy_degeneracies(t, PPS_Type(), use_result);
  }

  Sign
  check_easy_degeneracies(const Site_2& t, SSS_Type,
			  bool& use_result) const
  {
    CGAL_precondition( t.is_point() );

    use_result = false;

    // ADD THE CASES WHERE t IS AN ENDPOINT OF ONE OF THE SEGMENTS
    return ZERO;
  }

  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& st, PPP_Type) const
  {
    CGAL_precondition( st.is_point() );
    
    Point_2 t = st.point();

    Oriented_side os =
      side_of_oriented_circle(p_.point(), q_.point(), r_.point(), t);
    if ( os == ON_POSITIVE_SIDE ) { return NEGATIVE; }
    if ( os == ON_NEGATIVE_SIDE ) { return POSITIVE; }
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

  template<class Type>
  Sign incircle_p_no_easy(const Site_2& st, Type) const
  {
    CGAL_precondition( st.is_point() );

    FT r2 = squared_radius();

    Point_2 t = st.point();

    FT d2 = CGAL::square(x() - t.x()) +
      CGAL::square(y() - t.y());

    return CGAL::compare(d2, r2);
  }

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

    return CGAL::compare(d2, r2 * n2);
  }


  Site_2 other_site(const Site_2& p, const Site_2& seg) const
  {
    if (same_points(p, seg.source_site())){
      return  seg.target_site();
    }
    return  seg.source_site();
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& t, int) const
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

#ifndef CGAL_DISABLE_AM_CODE
    // code added by Andreas + Monique -- start
    if(v_type == PPP){
      Site_2 const *p1 = NULL;
      if(is_endpoint_of(p_, t)){
	p1 = &p_;
      } else if(is_endpoint_of(q_, t)){
	p1 = &q_;
      } else if(is_endpoint_of(r_, t)){
	p1 = &r_;
      }
      if(p1 != NULL){
        // As the Voronoi circle and the segment t touch in p1,
        // it is enough to check that the center and the non-touching point of the segment
        // are not in the same halfspace defined by the tangent line through p1
	Site_2 p2 = other_site(*p1, t);
	  Point_2 v(x(),y());
        
	  Compute_scalar_product_2 csp;
	  return -CGAL::sign( csp((v - p1->point()), (p2.point()- p1->point())) );
      }

    } else if(v_type == PPS){
      Site_2 const *p1, *p2, *seg;
      if(p_.is_point()){ 
	p1 = &p_;
	if(q_.is_point()){
	  p2 = &q_;
	  seg = &r_;
	} else {
	  p2 = &r_;
	  seg = &q_;
	} 
      } else {
	seg = &p_;
	p1 = &q_;
	p2 = &r_;
      }

      if(! is_endpoint_of(*p2, t)){
	std::swap(p1,p2);
      }
	if(is_endpoint_of(*p2, t)){
	  Site_2 tp = other_site(*p2, t);
	  Point_2 v(x(),y());
	  Compute_scalar_product_2 csp;
          return -CGAL::sign( csp((v - p2->point()), (tp.point()- p2->point())) );
	}
    }
    // code added by Andreas + Monique -- end
#endif // CGAL_DISABLE_AM_CODE

#ifndef CGAL_DISABLE_M_CODE
    // code added by Menelaos -- begin
    // in the code that follows we check whether one endpoint of the
    // query segment t is the same as the point p of a PSS circle. in
    // this case the result is known by taking the other point of t
    // and checking against the tangent to the Voronoi circle at p.
    if ( v_type == PSS ) {
      const Site_2* pt;
      if ( p_.is_point() ) { pt = &p_; }
      else if ( q_.is_point() ) { pt = &q_; }
      else { pt = &r_; }

      if ( is_endpoint_of(*pt, t) ) {
	Site_2 tp = other_site(*pt, t);
	Point_2 v(x(), y());
	Compute_scalar_product_2 csp;
	return
	  -CGAL::sign( csp(v - pt->point(), tp.point()- pt->point()) );
      }
    }
    // code added by Menelaos -- end
#endif // CGAL_DISABLE_M_CODE

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

    return incircle_s_no_easy(t, 0);
  }

  
  Sign incircle_s_no_easy(const Site_2& t, int) const
  {
    Sign d1, d2;
    if (  ( p_.is_point() && same_points(p_, t.source_site()) ) ||
	  ( q_.is_point() && same_points(q_, t.source_site()) ) ||
	  ( r_.is_point() && same_points(r_, t.source_site()) )  ) {
      d1 = ZERO;
    } else {
      d1 = incircle_p(t.source_site());
    }
    //if ( d1 == NEGATIVE ) { return NEGATIVE; }
    if ( certainly(d1 == NEGATIVE) ) { return NEGATIVE; }
    if (! is_certain(d1 == NEGATIVE) ) { return indeterminate<Sign>(); }

    if (  ( p_.is_point() && same_points(p_, t.target_site()) ) ||
	  ( q_.is_point() && same_points(q_, t.target_site()) ) ||
	  ( r_.is_point() && same_points(r_, t.target_site()) )  ) {
      d2 = ZERO;
    } else {
      d2 = incircle_p(t.target_site());
    }
    //if ( d2 == NEGATIVE ) { return NEGATIVE; }
    if (certainly( d2 == NEGATIVE ) ) { return NEGATIVE; }
    if (! is_certain( d2 == NEGATIVE ) ) { return indeterminate<Sign>(); }

    Line_2 l = compute_supporting_line(t.supporting_site());
    Sign sl = incircle(l);

    //if ( sl == POSITIVE ) { return sl; }
    if (certainly( sl == POSITIVE )) { return sl; }
    if (! is_certain( sl == POSITIVE )) { return indeterminate<Sign>(); }

    if ( sl == ZERO && (d1 == ZERO || d2 == ZERO) ) { return ZERO; }

    Oriented_side os1 = oriented_side(l, t.source());
    Oriented_side os2 = oriented_side(l, t.target());

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

    Sign s = incircle_s(t, 0);

    return s;
  }

  inline
  Sign incircle_s_no_easy(const Site_2& t) const
  {
    return incircle_s_no_easy(t, 0);
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
	Point_2 pref = p_ref().point();
	FT dx2 = CGAL::square(x() - pref.x());
	FT dy2 = CGAL::square(y() - pref.y());
	return dx2 + dy2;
      }
    case SSS:
      {
	Line_2 l = compute_supporting_line(p_.supporting_site());
	Homogeneous_point_2 q = compute_projection(l, point());

	FT dx2 = CGAL::square(x() - q.x());
	FT dy2 = CGAL::square(y() - q.y());
	return dx2 + dy2;
      }
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
  Voronoi_vertex_sqrt_field_C2(const Site_2& p,
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
      return incircle_p(t);
    }
    return incircle_s(t);
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
    return CGAL::sign(l.a() * x() + l.b() * y() + l.c());
  }

  Oriented_side oriented_side(const Line_2& l) const
  {
    return orientation(l);
  }

  //--------------------------------------------------------------------------

private:
  const Site_2& p_, q_, r_;

  vertex_t v_type;

  // index that indicates the refence point for the case PPS
  short pps_idx;

  FT ux, uy, uz;
};

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VEFTEX_SQRT_FIELD_C2_H

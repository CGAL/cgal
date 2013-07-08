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




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_SQRT_FIELD_NEW_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_SQRT_FIELD_NEW_C2_H



#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

template<class K>
class Voronoi_vertex_sqrt_field_new_C2
  : public Basic_predicates_C2<K>
{
public:
  typedef Basic_predicates_C2<K> Base;

  using Base::compute_supporting_line;
  using Base::oriented_side_of_line;
  using Base::opposite_line;
  using Base::compute_projection;

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
  // helpful methods
  //--------------------------------------------------------------------------

  // check whether the point p is an endpoint of the segment s
  inline
  bool is_endpoint_of(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );

    return ( same_points(p, s.source_site()) ||
	     same_points(p, s.target_site()) );
  }

  // given that p is an endpoint of seg (not checked), returns the
  // other endpoint of seg
  inline
  Site_2 other_site(const Site_2& sp, const Site_2& seg) const
  {
    CGAL_precondition( sp.is_point() && seg.is_segment() );

    if ( same_points(sp, seg.source_site()) ){
      return seg.target_site();
    }
    return seg.source_site();
  }



  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Voronoi vertex computation
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  // the Voronoi vertex of three points
  //--------------------------------------------------------------------------

  void
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
	     const PPP_Type&) const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
		       sr.is_point() );

    // the following check is not really needed in this
    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    Point_2 p = sp.point(), q = sq.point(), r = sr.point();

    FT np = CGAL::square(p.x()) + CGAL::square(p.y());
    FT nq = CGAL::square(q.x()) + CGAL::square(q.y());
    FT nr = CGAL::square(r.x()) + CGAL::square(r.y());

    FT ux, uy, uz;

    ux = np * (q.y() - r.y()) + nq * (r.y() - p.y()) + nr * (p.y() - q.y());
    uy = -(np * (q.x() - r.x()) + nq * (r.x() - p.x()) + nr * (p.x() - q.x()));
    uz = FT(2) * ( (q.x() * r.y() - r.x() * q.y()) +
		   (r.x() * p.y() - p.x() * r.y()) +
		   (p.x() * q.y() - q.x() * p.y()) );

    vv = Point_2(ux / uz, uy / uz);
  }


  //--------------------------------------------------------------------------
  // the Voronoi vertex of two points and a segment
  //--------------------------------------------------------------------------

  void
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
	     const PPS_Type&) const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
		       sr.is_segment() );

    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    FT a, b, c;
    compute_supporting_line(sr.supporting_site(), a, b, c);

    Point_2 pp = sp.point(), qq = sq.point();

    FT c_ = a * pp.x() + b * pp.y() + c;
    FT cq_ = a * qq.x() + b * qq.y() + c;

    if ( same_points(sp, sr.source_site()) ||
	 same_points(sp, sr.target_site()) ) {
      c_ = FT(0);
    }
    if ( same_points(sq, sr.source_site()) ||
	 same_points(sq, sr.target_site()) ) {
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

    Point_2 r_src = sr.source_site().point();
    Point_2 r_trg = sr.target_site().point();

    bool pq_xaligned = pp.y() == qq.y();
    bool r_xaligned = r_src.y() == r_trg.y();

    bool pq_yaligned = pp.x() == qq.x();
    bool r_yaligned = r_src.x() == r_trg.x();

    bool parallel = (pq_xaligned && r_xaligned) || (pq_yaligned && r_yaligned);
    // addition for parallel objects ends here

    if ( parallel || c_ == cq_ ) {
      FT e1 = CGAL::square(c_);
      FT J = nl * (a * n_ + FT(4) * c_ * x_) - FT(4) * a * e1;
      FT I = nl * (b * n_ + FT(4) * c_ * y_) - FT(4) * b * e1;
      FT X = FT(8) * nl * c_;

      // the homogeneous coordinates of the Voronoi vertex are:
      //          (J + pp.x() * X, I + pp.y() * X, X)

      vv = Point_2(J / X + pp.x(), I / X + pp.y());
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

    FT ux = J + pp.x() * X - FT(2) * y_ * sqrt_S;
    FT uy = I + pp.y() * X + FT(2) * x_ * sqrt_S;

    vv = Point_2(ux / X, uy / X);
  }


  //--------------------------------------------------------------------------
  // the Voronoi vertex of a point and two segments
  // also: the Voronoi vertex of a point and two lines
  //--------------------------------------------------------------------------


  void
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
	     const PSS_Type&) const
  {
#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      std::ofstream ofs("vv-failure-log.cin", std::ios_base::app);
      ofs.precision(16);
      ofs << sp << std::endl;
      ofs << sq << std::endl;
      ofs << sr << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

    CGAL_precondition( sp.is_point() && sq.is_segment() &&
		       sr.is_segment() );

    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    bool pq = is_endpoint_of(sp, sq);
    bool pr = is_endpoint_of(sp, sr);

    Point_2 pp = sp.point();

    if ( pq && pr ) {
      vv = pp;
      return;
    }


    FT a1, b1, c1, a2, b2, c2;
    compute_supporting_line(sq.supporting_site(), a1, b1, c1);
    compute_supporting_line(sr.supporting_site(), a2, b2, c2);

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
      if ( same_points(sp, sq.target_site()) ) {
	a1 = -a1;  b1 = -b1;  c1_ = -c1_;
      }
    }

    if ( sgn_c2_ == NEGATIVE ) {
      a2 = -a2;  b2 = -b2;  c2_ = -c2_;
    } else if ( sgn_c2_ == ZERO ) {

      CGAL_assertion( pr );
      if ( same_points(sp, sr.source_site()) ) {
	a2 = -a2;  b2 = -b2;  c2_ = -c2_;
      }
    }

    if ( pq ) {
      FT J = a1 * c2_;
      FT I = b1 * c2_;

      FT n1 = CGAL::square(a1) + CGAL::square(b1);
      FT n2 = CGAL::square(a2) + CGAL::square(b2);

      FT D1D2 = n1 * n2;

      FT uz = -a1 * a2 - b1 * b2 + CGAL::sqrt(D1D2);

      // the homogeneous coordinates of the Voronoi vertex are:
      //       (J + pp.x() * uz, uy = I + pp.y() * uz, uz)

      vv = Point_2(J / uz + pp.x(), I / uz + pp.y());

    } else if ( pr ) {
      FT J = a2 * c1_;
      FT I = b2 * c1_;

      FT n1 = CGAL::square(a1) + CGAL::square(b1);
      FT n2 = CGAL::square(a2) + CGAL::square(b2);

      FT D1D2 = n1 * n2;

      FT uz = -a1 * a2 - b1 * b2 + CGAL::sqrt(D1D2);

      // the homogeneous coordinates of the Voronoi vertex are:
      //      (J + pp.x() * uz, I + pp.y() * uz, uz)

      vv = Point_2(J / uz + pp.x(), I / uz + pp.y());

    } else {
      Line_2 lq(a1, b1, c1_);
      Line_2 lr(a2, b2, c2_);
      compute_pll(pp, lq, lr);
    }
  }


  void
  compute_pll(const Point_2& p, const Line_2& lq, const Line_2& lr) const
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

    FT uz = -a1a2 - b1b2 + CGAL::sqrt(D1D2);

    FT ux = J + p.x() * uz + sigma * CGAL::sqrt(u1);
    FT uy = I + p.y() * uz - rho * CGAL::sqrt(u2);

    vv = Point_2(ux / uz, uy / uz);
  }



  //--------------------------------------------------------------------------
  // the Voronoi vertex of three segments
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
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
	     const SSS_Type&) const
  {
    CGAL_precondition( sp.is_segment() && sq.is_segment() &&
		       sr.is_segment() );

    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    FT a[3], b[3], c[3];
    FT cx[3], cy[3], cz[3], sqrt_D[3];

    orient_lines(sp, sq, sr, a, b, c);

    for (int i = 0; i < 3; i++) {
      cx[i] = c[(i+1)%3] * b[(i+2)%3] - c[(i+2)%3] * b[(i+1)%3];
      cy[i] = -(c[(i+1)%3] * a[(i+2)%3] - c[(i+2)%3] * a[(i+1)%3]);
      cz[i] = -(a[(i+1)%3] * b[(i+2)%3] - a[(i+2)%3] * b[(i+1)%3]);

      FT d = CGAL::square(a[i]) + CGAL::square(b[i]);
      sqrt_D[i] = CGAL::sqrt(d);
    }

    FT ux = cx[0] * sqrt_D[0] + cx[1] * sqrt_D[1] + cx[2] * sqrt_D[2];
    FT uy = cy[0] * sqrt_D[0] + cy[1] * sqrt_D[1] + cy[2] * sqrt_D[2];
    FT uz = cz[0] * sqrt_D[0] + cz[1] * sqrt_D[1] + cz[2] * sqrt_D[2];

    vv = Point_2(ux / uz, uy / uz);
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Voronoi squared radius computation
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  template<class Type>
  inline
  FT
  squared_radius(const Point_2& vv,
		 const Site_2& p, const Site_2& /*q*/, const Site_2& /*r*/,
		 const Type&) const
  {
    CGAL_precondition( p.is_point() );

    Point_2 pp = p.point();
    FT dx2 = CGAL::square(vv.x() - pp.x());
    FT dy2 = CGAL::square(vv.y() - pp.y());

    return dx2 + dy2;
  }

  inline
  FT
  squared_radius(const Point_2& vv,
		 const Site_2& p, const Site_2& q, const Site_2& r,
		 const SSS_Type&) const
  {
    CGAL_USE(q);
    CGAL_USE(r);
    CGAL_assertion( p.is_segment() && q.is_segment() && r.is_segment() );

    Line_2 l = compute_supporting_line(p.supporting_site());
    Homogeneous_point_2 pref = compute_projection(l, vv);

    FT dx2 = CGAL::square(vv.x() - pref.x());
    FT dy2 = CGAL::square(vv.y() - pref.y());
    return dx2 + dy2;
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // the incircle test --- start
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  // given the Voronoi vertex vv or p, q and r, returns the result of
  // the incircle test when the query object t is a point
  template<class Type>
  Sign incircle_p(const Point_2& vv,
		  const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const Type& type) const
  {
    CGAL_precondition( t.is_point() );

    FT r2 = squared_radius(vv, p, q, r, type);
    Point_2 tt = t.point();
    FT d2 = CGAL::square(vv.x() - tt.x()) + CGAL::square(vv.y() - tt.y());

    return CGAL::compare(d2, r2);
  }

  //--------------------------------------------------------------------------
  // the first three objects are points and the query object is also a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const PPP_Type&) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_point() );
    CGAL_precondition( t.is_point() );
    
    Oriented_side os =
      side_of_oriented_circle(p.point(), q.point(), r.point(), t.point());
    if ( os == ON_POSITIVE_SIDE ) { return NEGATIVE; }
    if ( os == ON_NEGATIVE_SIDE ) { return POSITIVE; }
    return ZERO;
  }

  //--------------------------------------------------------------------------
  // the first three objects are two points and a segment and the query
  // object is a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const PPS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_point() );

    // easy degeneracies --- start

    // if t is one of p or q then we know the result which is ZERO
    if (  same_points(p, t) || same_points(q, t)  ) {
      return ZERO;
    }

    // if t is an endpoint of r, then t is necessarily outside the
    // Voronoi circle of p, q and r and thus the result is POSITIVE
    if ( is_endpoint_of(t, r) ) { return POSITIVE; }

    // easy degeneracies --- end

    compute_vv(p, q, r, type);
    return incircle_p(vv, p, q, r, t, type);
  }


  //--------------------------------------------------------------------------
  // the first three objects are a point and two segments and the query
  // object is a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const PSS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_point() );

    // easy degeneracies --- start

    // if p is a common point for q and r, then the Voronoi vertex of
    // p, q, and r is p. Since t cannot be the same as p, the result
    // is POSITIVE
    if ( is_endpoint_of(p, q) && is_endpoint_of(p, r) ) {
      return POSITIVE;
    }

    // if p and t are the same point, then return ZERO
    if ( same_points(p, t) ) { return ZERO; }
    
    // if t is an endpoint of either q or r, then t has to be outside
    // the Voronoi circle and thus the result is POSITIVE
    if ( is_endpoint_of(t, q) || is_endpoint_of(t, r) ) {
      return POSITIVE;
    }

    // easy degeneracies --- end

    compute_vv(p, q, r, type);
    return incircle_p(vv, p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  // the first three objects are segments and the query object is a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const SSS_Type& type) const
  {
    CGAL_precondition( p.is_segment() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_point() ); 

    // easy degeneracies --- start

    // if t is an endpoint of p, q or r, then t has to lie outside the
    // Voronoi circle of p, q and r and thus the result is positive
    if ( is_endpoint_of(t, p) || is_endpoint_of(t, q) ||
	 is_endpoint_of(t, r) ) {
      return POSITIVE;
    }

    // easy degeneracies --- end

    compute_vv(p, q, r, type);
    return incircle_p(vv, p, q, r, t, type);
  }


  //--------------------------------------------------------------------------
  // the incircle test when the query object is a line
  //--------------------------------------------------------------------------

  template<class Type>
  inline
  Sign
  incircle_xxxl(const Point_2& vv,
		const Site_2& p, const Site_2& q, const Site_2& r,
		const Line_2& l, const Type& type) const
  {
    FT r2 = squared_radius(vv, p, q, r, type);

    FT n2 = CGAL::square(l.a()) + CGAL::square(l.b());

    FT d2 = CGAL::square(l.a() * vv.x() + l.b() * vv.y() + l.c());

    return CGAL::compare(d2, r2 * n2);
  }


  inline
  Oriented_side
  oriented_side(const Point_2& vv, const Line_2& l, const Point_2& p) const
  {
    Line_2 l1(l.b(), -l.a(), l.a() * vv.y() - l.b() * vv.x());

    return oriented_side_of_line(l1, p);
  }


  //--------------------------------------------------------------------------
  // generic incircle test when the query object is a segment
  //--------------------------------------------------------------------------

  // first check is at least one of the endpoints of s is in conflict
  // with the Voronoi circle; in this case return NEGATIVE
  // otherwise test against the supporting line l of t; if l does not
  // conflict with the Voronoi circle return POSITIVE, otherwise check
  // if the endpoints of t are on the same oriented side of the line
  // perpendicular to l, passing through the Voronoi vertex of p, q,
  // and r, and respond accordingly: if they are on the same side
  // there is no conflict, otherwise there is a conflict.
  template<class Type>
  Sign incircle_xxxs(const Site_2& p, const Site_2& q, const Site_2& r,
		     const Site_2& t, const Type& type) const
  {
    CGAL_precondition( t.is_segment() );

    Sign d1, d2;
    if (  ( p.is_point() && same_points(p, t.source_site()) ) ||
	  ( q.is_point() && same_points(q, t.source_site()) ) ||
	  ( r.is_point() && same_points(r, t.source_site()) )  ) {
      d1 = ZERO;
    } else {
      d1 = incircle_p(p, q, r, t.source_site(), type);
    }

    if (  certainly(d1 == NEGATIVE)  ) { return NEGATIVE; }
    if (  !is_certain(d1 == NEGATIVE)  ) { return indeterminate<Sign>(); }

    if (   ( p.is_point() && same_points(p, t.target_site()) ) ||
	   ( q.is_point() && same_points(q, t.target_site()) ) ||
	   ( r.is_point() && same_points(r, t.target_site()) )  ) {
      d2 = ZERO;
    } else {
      d2 = incircle_p(p, q, r, t.target_site(), type);
    }

    if (  certainly( d2 == NEGATIVE )  ) { return NEGATIVE; }
    if (  !is_certain( d2 == NEGATIVE )  ) { return indeterminate<Sign>(); }

    Line_2 l = compute_supporting_line(t.supporting_site());
    compute_vv(p, q, r, type);
    Sign sl = incircle_xxxl(vv, p, q, r, l, type);
    
    if (  certainly( sl == POSITIVE )  ) { return sl; }
    if (  !is_certain( sl == POSITIVE )  ) { return indeterminate<Sign>(); }

    if ( sl == ZERO && (d1 == ZERO || d2 == ZERO) ) { return ZERO; }

    Oriented_side os1 = oriented_side(vv, l, t.source());
    Oriented_side os2 = oriented_side(vv, l, t.target());

    if ( sl == ZERO ) {
      if (os1 == ON_ORIENTED_BOUNDARY || os2 == ON_ORIENTED_BOUNDARY) {
	return ZERO;
      }
      return ( os1 == os2 ) ? POSITIVE : ZERO;
    }

    return (os1 == os2) ? POSITIVE : NEGATIVE;
  }


  //--------------------------------------------------------------------------
  // the first three objects are points and the query object is a segment
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const PPP_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_point() );
    CGAL_precondition( t.is_segment() );

    // easy degeneracies --- start

    // check if the endpoints of t are two of p, q and r
    unsigned int n_ends = 0;
    bool end_pt = is_endpoint_of(p, t);
    bool end_qt = is_endpoint_of(q, t);
    bool end_rt = is_endpoint_of(r, t);
    if ( end_pt ) ++n_ends;
    if ( end_qt ) ++n_ends;
    if ( end_rt ) ++n_ends;

    CGAL_assertion( n_ends < 3 );
    if ( n_ends == 2 ) { return NEGATIVE; }

#ifndef CGAL_DISABLE_AM_CODE
    // code added in previous version by Andreas + Monique -- start
    Site_2 const *pp1 = NULL;
    if ( end_pt ) pp1 = &p;
    else if ( end_qt ) pp1 = &q;
    else if ( end_rt ) pp1 = &r;
    if ( pp1 != NULL ) {
      // As the Voronoi circle and the segment t touch in p1,
      // it is enough to check that the center and the non-touching
      // point of the segment
      // are not in the same halfspace defined by the tangent line through p1
      Point_2 p1 = pp1->point();
      Point_2 p2 = other_site(*pp1, t).point();
      compute_vv(p, q, r, type);
        
      Compute_scalar_product_2 csp;
      return -CGAL::sign( csp(vv - p1, p2 - p1) );
    }
    // code added in previous version by Andreas + Monique -- end
#endif // CGAL_DISABLE_AM_CODE

    // easy degeneracies --- end

    return incircle_xxxs(p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  // the first three objects are two points and a segment and the
  // query object is a segment
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const PPS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_segment() );

    // easy degeneracies --- start

    // check if the endpoints of t are p and q
    bool end_pt = is_endpoint_of(p, t);
    bool end_qt = is_endpoint_of(q, t);

    if ( end_pt && end_qt ) { return NEGATIVE; }

#ifndef CGAL_DISABLE_AM_CODE
    // code added in previous version by Andreas + Monique -- start
    Site_2 const *pp1 = &p, *pp2 = &q;
    if ( !end_qt ) { std::swap(pp1, pp2); }

    if ( is_endpoint_of(*pp2, t) ) {
      Point_2 p1 = other_site(*pp2, t).point();
      Point_2 p2 = pp2->point();
      compute_vv(p, q, r, type);

      Compute_scalar_product_2 csp;
      return -CGAL::sign( csp(vv - p2, p1 - p2) );
    }
    // code added in previous version by Andreas + Monique -- end
#endif // CGAL_DISABLE_AM_CODE

    // easy degeneracies --- end

    return incircle_xxxs(p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  // the first three objects are a point and two segments and the query
  // object is a segment
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const PSS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_segment() );

    // easy degeneracies --- start

    // check if p is a common endpoint of q and r, in which case the
    // Voronoi circle degenerates to p
    if ( is_endpoint_of(p, q) && is_endpoint_of(p, r) ) {
      // case 1: the new segment is not adjacent to the center of the
      //         degenerate Voronoi circle, i.e., not adjacent to p
      if ( !is_endpoint_of(p, t) ) { return POSITIVE; }

      // check if t has the same support as either q or r
      if ( same_segments(q.supporting_site(), t.supporting_site()) ) {
	return ZERO;
      }

      if ( same_segments(r.supporting_site(), t.supporting_site()) ) {
	return ZERO;
      }


      Point_2 r_ = r.source(), q_ = q.source(), t_ = t.source();

      if ( same_points(q.source_site(), p) ) { q_ = q.target(); }
      if ( same_points(r.source_site(), p) ) { r_ = r.target(); }
      if ( same_points(t.source_site(), p) ) { t_ =  t.target(); }

      Point_2 p_ = p.point();

      if ( CGAL::orientation(p_, q_, t_) == LEFT_TURN &&
	   CGAL::orientation(p_, r_, t_) == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return ZERO;
    }

#ifndef CGAL_DISABLE_M_CODE
    // code added by Menelaos -- begin

    // in the code that follows we check whether one endpoint of the
    // query segment t is the same as the point p of a PSS circle. in
    // this case the result is known by taking the other point of t
    // and checking against the tangent to the Voronoi circle at p.
    if ( is_endpoint_of(p, t) ) {
      Point_2 p1 = p.point();
      Point_2 p2 = other_site(p, t).point();
      compute_vv(p, q, r, type);
      Compute_scalar_product_2 csp;
      return -CGAL::sign( csp(vv - p1, p2 - p1) );
    }
    // code added by Menelaos -- end
#endif // CGAL_DISABLE_M_CODE

    // check if t has the same support as either q or r
    if ( same_segments(q.supporting_site(), t.supporting_site()) ) {
      return POSITIVE;
    }

    if ( same_segments(r.supporting_site(), t.supporting_site()) ) {
      return POSITIVE;
    }

    // easy degeneracies --- end

    return incircle_xxxs(p, q, r, t, type);
  }


  //--------------------------------------------------------------------------
  // the first three objects are segments and the query object is also
  // a segment
  //--------------------------------------------------------------------------

  inline
  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t, const SSS_Type& type) const
  {
    CGAL_precondition( p.is_segment() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_segment() );

    return incircle_xxxs(p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // the incircle test --- end
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------



  vertex_t
  compute_type(const Site_2& s1, const Site_2& s2, const Site_2& s3) const
  {
    int npts = 0;
    if ( s1.is_point() ) ++npts;
    if ( s2.is_point() ) ++npts;
    if ( s3.is_point() ) ++npts;

    switch ( npts ) {
    case 0:
      return SSS;
      break;
    case 1:
      return PSS;
      break;
    case 2:
      return PPS;
      break;
    default:
      return PPP;
    }
  }


  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_point() );

    switch ( v_type ) {
    case PPP:
      return incircle_p(p, q, r, t, PPP_Type());
    case PPS:
      if ( p.is_segment() ) {
	return incircle_p(q, r, p, t, PPS_Type());
      } else if ( q.is_segment() ) {
	return incircle_p(r, p, q, t, PPS_Type());
      } else {
	return incircle_p(p, q, r, t, PPS_Type());
      }
    case PSS:
      if ( p.is_point() ) {
	return incircle_p(p, q, r, t, PSS_Type());
      } else if ( q.is_point() ) {
	return incircle_p(q, r, p, t, PSS_Type());
      } else {
	return incircle_p(r, p, q, t, PSS_Type());
      }
      return incircle_p(p, q, r, t, PSS_Type());
    default: // case SSS
      return incircle_p(p, q, r, t, SSS_Type());
    }
  }



  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_segment() );

    switch ( v_type ) {
    case PPP:
      return incircle_s(p, q, r, t, PPP_Type());
    case PPS:
      if ( p.is_segment() ) {
	return incircle_s(q, r, p, t, PPS_Type());
      } else if ( q_.is_segment() ) {
	return incircle_s(r, p, q, t, PPS_Type());
      } else {
	return incircle_s(p, q, r, t, PPS_Type());
      }
    case PSS:
      if ( p.is_point() ) {
	return incircle_s(p, q, r, t, PSS_Type());
      } else if ( q.is_point() ) {
	return incircle_s(q, r, p, t, PSS_Type());
      } else {
	return incircle_s(r, p, q, t, PSS_Type());
      }
    default: // case SSS
      return incircle_s(p, q, r, t, SSS_Type());
    }
  }


public:
  Voronoi_vertex_sqrt_field_new_C2(const Site_2& p,
				   const Site_2& q,
				   const Site_2& r)
    : p_(p), q_(q), r_(r), is_vv_computed(false)
  {
    v_type = compute_type(p, q, r);
  }


  inline bool is_degenerate_Voronoi_circle() const
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


  Point_2 degenerate_point() const
  {
    CGAL_precondition( is_degenerate_Voronoi_circle() );
    if ( p_.is_point() ) return p_.point();
    if ( q_.is_point() ) return q_.point();
    return r_.point();
  }


  Point_2 point() const
  {
    if ( is_degenerate_Voronoi_circle() ) {
      return degenerate_point();
    }

    if ( !is_vv_computed ) {
      switch ( v_type ) {
      case PPP:
	compute_vv(p_, q_, r_, PPP_Type());
	break;
      case PPS:
	if ( p_.is_segment() ) {
	  compute_vv(q_, r_, p_, PPS_Type());
	} else if ( q_.is_segment() ) {
	  compute_vv(r_, p_, q_, PPS_Type());
	} else {
	  compute_vv(p_, q_, r_, PPS_Type());
	}
	break;
      case PSS:
	if ( p_.is_point() ) {
	  compute_vv(p_, q_, r_, PSS_Type());
	} else if ( q_.is_point() ) {
	  compute_vv(q_, r_, p_, PSS_Type());
	} else {
	  compute_vv(r_, p_, q_, PSS_Type());
	}
	break;
      default: // case SSS:
	compute_vv(p_, q_, r_, SSS_Type());
	break;
      }
    }

    return vv;
  }


  inline Sign incircle(const Site_2& t) const
  {
    if ( t.is_point() ) {
      return incircle_p(p_, q_, r_, t);
    }
    return incircle_s(p_, q_, r_, t);
  }


  inline Sign operator()(const Site_2& p, const Site_2& q, const Site_2& r,
			 const Site_2& t) const
  {
    if ( t.is_point() ) {
      return incircle_p(p, q, r, t);
    }
    return incircle_s(p, q, r, t);
  }

private:

  template<class Type>
  Sign incircle_p_no_easy(const Point_2& vv,
			  const Site_2& p, const Site_2& q, const Site_2& r,
			  const Site_2& t, const Type& type) const
  {
    CGAL_precondition( t.is_point() );

    FT r2 = squared_radius(vv, p, q, r, type);

    Point_2 tt = t.point();

    FT d2 = CGAL::square(vv.x() - tt.x()) + CGAL::square(vv.y() - tt.y());

    return CGAL::compare(d2, r2);
  }

  
  Sign incircle_p_no_easy(const Site_2& p, const Site_2& q, const Site_2& r,
			  const Site_2& t) const
  {
    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_p(p, q, r, t, PPP_Type());
      break;
    case PPS:
      PPS_Type pps;
      if ( p.is_segment() ) {
	compute_vv(q, r, p, pps);
	s = incircle_p_no_easy(vv, q, r, p, t, pps);
      } else if ( q.is_segment() ) {
	compute_vv(r, p, q, pps);
	s = incircle_p_no_easy(vv, r, p, q, t, pps);
      } else {
	compute_vv(p, q, r, pps);
	s = incircle_p_no_easy(vv, p, q, r, t, pps);
      }
      break;
    case PSS:
      PSS_Type pss;
      if ( p.is_point() ) {
	compute_vv(p, q, r, pss);
	s = incircle_p_no_easy(vv, p, q, r, t, pss);
      } else if ( q.is_point() ) {
	compute_vv(q, r, p, pss);
	s = incircle_p_no_easy(vv, q, r, p, t, pss);
      } else {
	compute_vv(r, p, q, pss);
	s = incircle_p_no_easy(vv, r, p, q, t, pss);
      }
      break;
    case SSS:
      SSS_Type sss;
      compute_vv(p, q, r, sss);
      s = incircle_p_no_easy(vv, p, q, r, t, sss);
      break;
    }

    return s;
  }


  Sign incircle_s_no_easy(const Site_2& p, const Site_2& q, const Site_2& r,
			  const Site_2& t) const 
  {
    switch ( v_type ) {
    case PPP:
      return incircle_xxxs(p, q, r, t, PPP_Type());
    case PPS:
      if ( p.is_segment() ) {
	return incircle_xxxs(q, r, p, t, PPS_Type());
      } else if ( q_.is_segment() ) {
	return incircle_xxxs(r, p, q, t, PPS_Type());
      } else {
	return incircle_xxxs(p, q, r, t, PPS_Type());
      }
    case PSS:
      if ( p.is_point() ) {
	return incircle_xxxs(p, q, r, t, PSS_Type());
      } else if ( q.is_point() ) {
	return incircle_xxxs(q, r, p, t, PSS_Type());
      } else {
	return incircle_xxxs(r, p, q, t, PSS_Type());
      }
    default: // case SSS:
      return incircle_xxxs(p, q, r, t, SSS_Type());
    }
  }




public:
  inline Sign incircle_no_easy(const Site_2& t) const
  {
    Sign s;

    if ( t.is_point() ) {
      s = incircle_p_no_easy(p_, q_, r_, t);
    } else {
      s = incircle_s_no_easy(p_, q_, r_, t);
    }

    return s;
  }


  Orientation orientation(const Line_2& l) const
  {
    switch ( v_type ) {
    case PPP:
      compute_vv(p_, q_, r_, PPP_Type());
      break;
    case PPS:
      if ( p_.is_segment() ) {
	compute_vv(q_, r_, p_, PPS_Type());
      } else if ( q_.is_segment() ) {
	compute_vv(r_, p_, q_, PPS_Type());
      } else {
	compute_vv(p_, q_, r_, PPS_Type());
      }
      break;
    case PSS:
      if ( p_.is_point() ) {
	compute_vv(p_, q_, r_, PSS_Type());
      } else if ( q_.is_point() ) {
	compute_vv(q_, r_, p_, PSS_Type());
      } else {
	compute_vv(r_, p_, q_, PSS_Type());
      }
      break;
    default: // case SSS:
      compute_vv(p_, q_, r_, SSS_Type());
      break;
    }

    return CGAL::sign(l.a() * vv.x() + l.b() * vv.y() + l.c());
  }


  inline Oriented_side oriented_side(const Line_2& l) const
  {
    return orientation(l);
  }

private:
  // the defining sites of the Voronoi vertex
  const Site_2& p_, &q_, &r_;

  // indicates whether the Voronoi vertex has been computed
  mutable bool is_vv_computed; 

  // the type of the Voronoi vertex
  vertex_t v_type;

  // the computed Voronoi vertex is cached in this variable
  mutable Point_2 vv;
};





} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_VORONOI_VERTEX_SQRT_FIELD_NEW_C2_H

// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/predicates/Apollonius_graph_predicates_C2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_PREDICATES_C2_H
#define CGAL_APOLLONIUS_GRAPH_PREDICATES_C2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>
#include <CGAL/Point_2.h>
#include <CGAL/Weighted_point.h>

#include <CGAL/predicates/Apollonius_graph_predicate_constructions_C2.h>

#include <CGAL/functions_on_signs.h>
#include <CGAL/predicates/compare_quadratic.h>

CGAL_BEGIN_NAMESPACE


//--------------------------------------------------------------------

template<class Point, class Weight>
class Sign_of_Voronoi_radius
{
public:
  typedef Voronoi_radius<Point, Weight>   Voronoi_radius;
  typedef typename Point::R::FT           FT;

public:

  inline
  Sign operator()(const Voronoi_radius& vr, Sqrt_field_tag )
    {
      FT r = vr.c1() + vr.c2() * CGAL_NTS sqrt(vr.delta());
      return CGAL_NTS sign(r);
    }


  inline
  Sign operator()(const Voronoi_radius& vr, Ring_tag )
    {
      // this is another way of doing this; the degree becomes 10
      // instead of 5 in this case.
      //      return sign_a_plus_b_x_sqrt_c(vr.c1(), vr.c2(), vr.delta());

      bool is_first_root = CGAL_NTS is_negative(vr.c2());

      Sign s_beta = CGAL_NTS sign(vr.b());
      Sign s_gamma = CGAL_NTS sign(vr.c());

      // the existence test
      Sign sign_of_root;
      if ( is_first_root ) {
	sign_of_root = sign_of_first_root(s_beta, s_gamma);
      } else {
	sign_of_root = sign_of_second_root(s_beta, s_gamma);
      }
      return sign_of_root;
    }
};

template<class Point, class Weight>
class Bounded_side_of_CCW_circle
{
private:
  Sign_of_Voronoi_radius<Point, Weight> test;
public:
  typedef Voronoi_radius<Point, Weight>   Voronoi_radius;
  typedef typename Point::R::FT           FT;
public:

  template<class Tag>
  Bounded_side
  operator()(const Voronoi_radius& vr, Tag tag)
    {
      return Bounded_side(  opposite( test(vr, tag ) )  );
    }
};


//--------------------------------------------------------------------

template<class Point, class Weight>
class Sign_of_distance_from_bitangent_line
{
public:
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef typename Point::R::FT                   FT;

public:

  Sign
  operator()(const Bitangent_line& bl, const Weighted_point& q,
	     Sqrt_field_tag)
    {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::distance_from_bitangent_counter++;
#endif
      FT a = bl.a1() + bl.a2() * CGAL_NTS sqrt(bl.delta());
      FT b = bl.b1() + bl.b2() * CGAL_NTS sqrt(bl.delta());
      FT c = bl.c1() + bl.c2() * CGAL_NTS sqrt(bl.delta());
      FT r = a * q.x() + b * q.y() + c - q.weight() * bl.d();
      return CGAL_NTS sign(r);
    }

  Sign
  operator()(const Bitangent_line& bl, const Weighted_point& q, Ring_tag)
    {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::distance_from_bitangent_counter++;
#endif
      FT A = bl.a1() * q.x() + bl.b1() * q.y() + bl.c1()
	- q.weight() * bl.d();
      FT B = bl.a2() * q.x() + bl.b2() * q.y() + bl.c2();
      return sign_a_plus_b_x_sqrt_c(A, B, bl.delta());
    }
};

//--------------------------------------------------------------------


template<class Point, class Weight>
class Sign_of_distance_from_CCW_circle
{
public:
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef typename Point::R::FT                   FT;
public:

  Sign
  operator()(const Bitangent_line& bl,
	     const Inverted_weighted_point& v, Sqrt_field_tag)
    {
      FT a = bl.a1() + bl.a2() * CGAL_NTS sqrt(bl.delta());
      FT b = bl.b1() + bl.b2() * CGAL_NTS sqrt(bl.delta());
      FT c = bl.c1() + bl.c2() * CGAL_NTS sqrt(bl.delta());
      FT r = a * v.x() + b * v.y() + c * v.p() - v.weight() * bl.d();
      return CGAL_NTS sign(r);
    }

  Sign
  operator()(const Bitangent_line& bl,
	     const Inverted_weighted_point& v, Ring_tag)
    {
      FT A = bl.a1() * v.x() + bl.b1() * v.y() + bl.c1() * v.p()
	- v.weight() * bl.d();
      FT B = bl.a2() * v.x() + bl.b2() * v.y() + bl.c2() * v.p();

      return sign_a_plus_b_x_sqrt_c(A, B, bl.delta());
    }
};


template < class Weighted_point >
class Weighted_point_less_than
{
public:
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2) const
  {
    if ( p1.x() == p2.x() ) {
      return p1.y() < p2.y();
    }
    return p1.x() < p2.x();
  }
};

template < class R, class W = typename R::RT >
class Incircle_test
{
public:
  typedef Point_2<R>                              Point;
  typedef W                                       Weight;

  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Weighted_point_inverter<Point, Weight>  Weighted_point_inverter;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef Voronoi_radius<Point, Weight>           Voronoi_radius;
  typedef typename R::FT                          FT;

  typedef Bounded_side_of_CCW_circle<Point, Weight>
                                               Bounded_side_of_CCW_circle;
  typedef Sign_of_distance_from_bitangent_line<Point, Weight>
                                     Sign_of_distance_from_bitangent_line;
  typedef Sign_of_distance_from_CCW_circle<Point, Weight>
                                         Sign_of_distance_from_CCW_circle;

private:
  Orientation
  orientation(const Bitangent_line& l, const Point& p,
	      Sqrt_field_tag) const
    {
      FT A = l.a1() * p.x() + l.b1() * p.y() + l.c1();
      FT B = l.a2() * p.x() + l.b2() * p.y() + l.c2();
      FT P = A + B * CGAL_NTS sqrt(l.delta());
      return CGAL_NTS sign(P);
    }

  Orientation
  orientation(const Bitangent_line& l, const Point& p,
	      Ring_tag) const
    {
      FT A = l.a1() * p.x() + l.b1() * p.y() + l.c1();
      FT B = l.a2() * p.x() + l.b2() * p.y() + l.c2();
      return sign_a_plus_b_x_sqrt_c(A, B, l.delta());
    }

  
  Orientation
  orientation(const Bitangent_line& l,
	      const Inverted_weighted_point& u) const
    {
      FT A = l.a1() * u.x() / u.p() + l.b1() * u.y() / u.p() + l.c1();
      FT B = l.a2() * u.x() / u.p() + l.b2() * u.y() / u.p() + l.c2();
      FT P = A + B * CGAL_NTS sqrt(l.delta());
      return CGAL_NTS sign(P);
    }

public:

  template<class Method_tag>
  inline Sign operator()(const Weighted_point& p1,
			 const Weighted_point& p2,
			 const Weighted_point& p3,
			 const Weighted_point& q,
			 Method_tag tag) const {
#ifdef AG2_PROFILE_PREDICATES
    ag2_predicate_profiler::incircle_counter++;
#endif
    //
    Weighted_point_inverter inverter(p1);
    Inverted_weighted_point u2 = inverter(p2);
    Inverted_weighted_point u3 = inverter(p3);

    Voronoi_radius vr_123(u2, u3);

    Bounded_side bs = Bounded_side_of_CCW_circle()(vr_123, tag );

    if ( bs != ON_UNBOUNDED_SIDE ) { return NEGATIVE; }

    Inverted_weighted_point v = inverter(q);
    Bitangent_line blinv_23(u2, u3);
    Sign s = Sign_of_distance_from_CCW_circle()(blinv_23, v, tag);
    return s;
  }

  template<class Method_tag>
  inline Sign operator()(const Weighted_point& p1,
			 const Weighted_point& p2,
			 const Weighted_point& q,
			 Method_tag tag) const {
    //
    Bitangent_line bl_21(p2, p1);
    Sign s = Sign_of_distance_from_bitangent_line()(bl_21, q, tag);
    if ( s != ZERO ) { return s; }

    Bitangent_line bl1_perp = bl_21.perpendicular(p1.point());
    Bitangent_line bl2_perp = bl_21.perpendicular(p2.point());
    Orientation o1 = orientation(bl1_perp, q.point(), tag);
    Orientation o2 = orientation(bl2_perp, q.point(), tag);

    CGAL_assertion( o1 != COLLINEAR && o2 != COLLINEAR );
    if ( o1 == o2 ) { return POSITIVE; }
    return NEGATIVE;
  }
};

//--------------------------------------------------------------------

template<class Point, class Weight>
class Orientation_wrt_symmetry_axis
{
public:
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef typename Point::R::FT                   FT;

public:

  Orientation
  operator()(const Voronoi_circle& vc, const Point& p1,
	     const Point& p2, Sqrt_field_tag)
    {
      FT a = vc.a1() + vc.a2() * CGAL_NTS sqrt(vc.delta());
      FT b = vc.b1() + vc.b2() * CGAL_NTS sqrt(vc.delta());
      FT det = a * (p2.y() - p1.y()) - b * (p2.x() - p1.x());
      return CGAL_NTS sign(det);
    }

  Orientation
  operator()(const Voronoi_circle& vc, const Point& p1,
	     const Point& p2, Ring_tag)
    {
      FT dx = p2.x() - p1.x();
      FT dy = p2.y() - p1.y();
      FT A = vc.a1() * dy - vc.b1() * dx;
      FT B = vc.a2() * dy - vc.b2() * dx;
      return sign_a_plus_b_x_sqrt_c(A, B, vc.delta());
    }
};


//--------------------------------------------------------------------

template<class Point, class Weight>
class Compare_Voronoi_radii
{
public:
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef typename Point::R::FT                   FT;

private:

  Sign sign_of_P4(const FT& u, const FT& v,
		  const FT& du, const FT& dv, const FT& dr,
		  const FT& Du, const FT& Dv, const FT& Dr)
    {
      std::pair<FT,FT> factors = factors_of_P4(u, v, du, dv, dr,
					       Du, Dv, Dr);
      Sign s1 = CGAL_NTS sign(factors.first);
      Sign s2 = CGAL_NTS sign(factors.second);
      return (s1 * s2);
    }

  std::pair<FT,FT>
  factors_of_P4(const FT& u, const FT& v,
		const FT& du, const FT& dv, const FT& dr,
		const FT& Du, const FT& Dv, const FT& Dr)
    {
      FT u2 = CGAL_NTS square(u);
      FT v2 = CGAL_NTS square(v);

      FT du2 = CGAL_NTS square(du);
      FT Du2 = CGAL_NTS square(Du);

      FT dv2 = CGAL_NTS square(dv);
      FT Dv2 = CGAL_NTS square(Dv);

      FT dr2 = CGAL_NTS square(dr);
      FT Dr2 = CGAL_NTS square(Dr);

      FT u2_P_v2 = u2 + v2;
      FT u2_M_v2 = u2 - v2;
      FT uv = FT(2) * u * v;

      FT drDr = FT(2) * dr * Dr;

      FT du2_P_dv2 = du2 + dv2;
      FT Du2_P_Dv2 = Du2 + Dv2;

      FT uU_P_vV = du * Du + dv * Dv;
      FT uU_M_vV = du * Du - dv * Dv;

      FT uV_P_Uv = du * Dv + Du * dv;
      FT uV_M_Uv = du * Dv - Du * dv;
      

      FT F1 = du2_P_dv2 * Dr2 + Du2_P_Dv2 * dr2
	- uU_P_vV * drDr - CGAL_NTS square(uV_M_Uv);

      FT F2 = CGAL_NTS square(u2_P_v2) * (du2_P_dv2 * Dr2 + Du2_P_Dv2 * dr2);
      F2 -= u2_P_v2 * (u2_M_v2 * uU_M_vV + uv * uV_P_Uv) * drDr;
      F2 -= CGAL_NTS square(u2_M_v2 * uV_P_Uv - uv * uU_M_vV);

      std::pair<FT, FT> factors(F1,F2);
      return factors;
    }

public:
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     Sqrt_field_tag)
    {
      FT c1 = (vc1.c1() + vc1.c2() * CGAL_NTS sqrt(vc1.delta())) / vc1.d();
      FT c2 = (vc2.c1() + vc2.c2() * CGAL_NTS sqrt(vc2.delta())) / vc2.d();

      Comparison_result r = CGAL_NTS compare(c2, c1);
      return r;
    }

  // this is the naive way but without divisions and square roots; the
  // degree becomes 36 in this case.
  /*
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     Ring_tag)
    {
      FT A = vc1.c1() * vc2.d() - vc2.c1() * vc1.d();
      FT B = vc1.c2() * vc2.d();
      FT C = -vc2.c2() * vc1.d();
      FT E = vc1.delta();
      FT F = vc2.delta();

      Sign s = sign_a_plus_b_x_sqrt_e_plus_c_x_sqrt_f(A,B,C,E,F);

      if ( s == ZERO ) { return EQUAL; }
      return ( s == POSITIVE ) ? SMALLER : LARGER;
    }
  */

  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     Ring_tag)
    {
      bool is_first_root1 = vc1.is_first_root();
      bool is_first_root2 = vc2.is_first_root();

      CGAL_precondition( CGAL_NTS is_positive(vc1.alpha()) );
      CGAL_precondition( CGAL_NTS is_positive(vc2.alpha()) );

      Comparison_result r;
      if ( is_first_root1 && is_first_root2 ) {
	r = ke_compare_l1_l2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			     vc2.alpha(), vc2.beta(), vc2.gamma());
      } else if ( is_first_root1 && !is_first_root2 ) {
	r = ke_compare_l1_r2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			     vc2.alpha(), vc2.beta(), vc2.gamma());
      } else if ( !is_first_root1 && is_first_root2 ) {
	r = ke_compare_r1_l2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			     vc2.alpha(), vc2.beta(), vc2.gamma());
      } else {
	r = ke_compare_r1_r2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			     vc2.alpha(), vc2.beta(), vc2.gamma());
      }

#ifdef COMPARATOR_PROFILER
      if ( comparator_profiler::count_cases ) {
	// count cases only for the tree r1-r2
	if ( !is_first_root1 && !is_first_root2 ) {
	  comparator_profiler::count_case(vc1.alpha(), vc1.beta(),
					  vc1.gamma(),
					  vc2.alpha(), vc2.beta(),
					  vc2.gamma());
	}
      }
#endif

      if ( r == EQUAL ) { return EQUAL; }
      return ( r == LARGER ) ? SMALLER : LARGER;
    }

  // this uses the DFMT trees; slightly slower but same degree (20).
  /*
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     Ring_tag)
    {
      bool is_first_root1 = vc1.is_first_root();
      bool is_first_root2 = vc2.is_first_root();

      CGAL_precondition( CGAL_NTS is_positive(vc1.alpha()) );
      CGAL_precondition( CGAL_NTS is_positive(vc2.alpha()) );

      Comparison_result r;
      if ( is_first_root1 && is_first_root2 ) {
	r = dfmt_compare_l1_l2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			  vc2.alpha(), vc2.beta(), vc2.gamma());
      } else if ( is_first_root1 && !is_first_root2 ) {
	r = dfmt_compare_l1_r2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			  vc2.alpha(), vc2.beta(), vc2.gamma());
      } else if ( !is_first_root1 && is_first_root2 ) {
	r = dfmt_compare_r1_l2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			  vc2.alpha(), vc2.beta(), vc2.gamma());
      } else {
	r = dfmt_compare_r1_r2(vc1.alpha(), vc1.beta(), vc1.gamma(),
			  vc2.alpha(), vc2.beta(), vc2.gamma());
      }

      if ( r == EQUAL ) { return EQUAL; }
      return ( r == LARGER ) ? SMALLER : LARGER;
    }
  */
};


//--------------------------------------------------------------------

template<class Point, class Weight>
class Order_on_finite_bisector
{
public:
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef typename Point::R::FT                   FT;
  typedef Orientation_wrt_symmetry_axis<Point, Weight>
                                   Orientation_wrt_symmetry_axis;
  typedef Compare_Voronoi_radii<Point, Weight>
                                           Compare_Voronoi_radii;
  
public:
  template<class Method_tag>
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     const Weighted_point& p1, const Weighted_point& p2,
	     Method_tag tag)
    {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::order_on_bisector_counter++;
#endif

      Orientation o1 =
	Orientation_wrt_symmetry_axis()(vc1, p1, p2, tag);
      Orientation o2 =
	Orientation_wrt_symmetry_axis()(vc2, p1, p2, tag);

      Comparison_result cr;
      if ( o1 == LEFT_TURN ) {
	if ( o2 != LEFT_TURN ) { return SMALLER; }
	Comparison_result r = Compare_Voronoi_radii()(vc1, vc2, tag);

	if ( r == EQUAL ) {
	  cr = EQUAL;
	} else {
	  cr = (r == LARGER ) ? SMALLER : LARGER;
	}
      } else if ( o1 == COLLINEAR ) {
	if ( o2 == COLLINEAR ) {
	  cr = EQUAL;
	} else {
	  cr = (o2 == LEFT_TURN) ? LARGER : SMALLER;
	}
      } else {
	if ( o2 != RIGHT_TURN ) {
	  cr = LARGER;
	} else {
	  Comparison_result r =
	    Compare_Voronoi_radii()(vc1, vc2, tag);
	  cr = r;
	}
      }

      return cr;
    }
};


//--------------------------------------------------------------------

template < class R, class W = typename R::RT >
class Finite_edge_test
{
public:
  typedef Point_2<R>                              Point;
  typedef W                                       Weight;

  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Weighted_point_inverter<Point, Weight>  Weighted_point_inverter;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef Voronoi_radius<Point, Weight>           Voronoi_radius;
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef typename R::FT                          FT;

  typedef Bounded_side_of_CCW_circle<Point, Weight>
                                               Bounded_side_of_CCW_circle;
  typedef Sign_of_distance_from_bitangent_line<Point, Weight>
                                     Sign_of_distance_from_bitangent_line;
  typedef Sign_of_distance_from_CCW_circle<Point, Weight>
                                         Sign_of_distance_from_CCW_circle;
  typedef Order_on_finite_bisector<Point, Weight>
                                                 Order_on_finite_bisector;

public:
  template<class Method_tag>
  bool
  operator()(const Weighted_point& p1,
	     const Weighted_point& p2,
	     const Weighted_point& p3,
	     const Weighted_point& p4,
	     const Weighted_point& q, bool b, Method_tag tag) {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::shadow_region_type_counter++;
#endif
    //
    Weighted_point_inverter inverter(p1);
    Inverted_weighted_point u2 = inverter(p2);
    Inverted_weighted_point v = inverter(q);
    //
    Voronoi_radius vr_12q(u2, v);
    Voronoi_radius vr_1q2 = vr_12q.get_symmetric();

    Bounded_side bs1 = Bounded_side_of_CCW_circle()(vr_12q, tag );
    Bounded_side bs2 = Bounded_side_of_CCW_circle()(vr_1q2, tag );

    bool is_bs1 = (bs1 == ON_UNBOUNDED_SIDE);
    bool is_bs2 = (bs2 == ON_UNBOUNDED_SIDE);

    // both the ccw and cw circles do not exist
    if ( !is_bs1 && !is_bs2 ) {
      return b;
    }

    // the ccw circle exists but not the cw
    if ( is_bs1 && !is_bs2 ) {
      return b;
    }

    // the cw circle exists but not the ccw
    if ( !is_bs1 && is_bs2 ) {
      return b;
    }


    // both circles exist

    // check whether the shadow region is connected, i.e., wether it is
    // of the form (a, b) or (-oo, a) U (b, +oo)
    Bitangent_line bl_12(p1, p2);

    Sign stc =
      Sign_of_distance_from_bitangent_line()(bl_12, q, tag);

    CGAL_assertion( stc != ZERO );
    bool is_shadow_region_connected = (stc == POSITIVE);

    if ( is_shadow_region_connected ) {
      if ( b ) { return true; }

      Inverted_weighted_point u3 = inverter(p3);
      Bitangent_line blinv_23(u2, u3);

      Voronoi_circle vc_123(blinv_23);
      Voronoi_circle vc_12q(vr_12q);


      Comparison_result r =
	Order_on_finite_bisector()(vc_123, vc_12q, p1, p2, tag);

      if ( r != SMALLER ) { return false; }

      Inverted_weighted_point u4 = inverter(p4);
      Bitangent_line blinv_42(u4, u2);

      Voronoi_circle vc_142(blinv_42);
      Voronoi_circle vc_1q2(vr_1q2);
      r = Order_on_finite_bisector()(vc_142, vc_1q2, p1, p2, tag);

      return ( r == LARGER );
    }

    // the shadow region is of the form (-oo, a) U (b, +oo)
    if ( !b ) { return false; }

    Inverted_weighted_point u3 = inverter(p3);
    Bitangent_line blinv_23(u2, u3);

    Voronoi_circle vc_123(blinv_23);
    Voronoi_circle vc_1q2(vr_1q2);
    Comparison_result r =
      Order_on_finite_bisector()(vc_123, vc_1q2, p1, p2, tag);

    if ( r != SMALLER ) { return true; }

    Inverted_weighted_point u4 = inverter(p4);
    Bitangent_line blinv_42(u4, u2);

    Voronoi_circle vc_142(blinv_42);
    Voronoi_circle vc_12q(vr_12q);
    r = Order_on_finite_bisector()(vc_142, vc_12q, p1, p2, tag);

    return ( r != LARGER );
  }
};

//--------------------------------------------------------------------

template < class R, class W = typename R::RT >
class Finite_edge_test_degenerated
{
public:
  typedef Point_2<R>                              Point;
  typedef W                                       Weight;

  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Weighted_point_inverter<Point, Weight>  Weighted_point_inverter;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef Voronoi_radius<Point, Weight>           Voronoi_radius;
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef typename R::FT                          FT;

  typedef Bounded_side_of_CCW_circle<Point, Weight>
                                               Bounded_side_of_CCW_circle;
  typedef Sign_of_distance_from_bitangent_line<Point, Weight>
                                     Sign_of_distance_from_bitangent_line;
  typedef Sign_of_distance_from_CCW_circle<Point, Weight>
                                         Sign_of_distance_from_CCW_circle;
  typedef Order_on_finite_bisector<Point, Weight>
                                                 Order_on_finite_bisector;
public:

  template<class Method_tag>
  bool
  operator()(const Weighted_point& p1,
	     const Weighted_point& p2,
	     const Weighted_point& p3,
	     const Weighted_point& q, bool b, Method_tag tag) {
#ifdef AG2_PROFILE_PREDICATES
    ag2_predicate_profiler::shadow_region_type_counter++;
#endif
    //
    Weighted_point_inverter inverter(p1);
    Inverted_weighted_point u2 = inverter(p2);
    Inverted_weighted_point v = inverter(q);

    Voronoi_radius vr_12q(u2, v);
    Voronoi_radius vr_1q2 = vr_12q.get_symmetric();

    Bounded_side bs1 = Bounded_side_of_CCW_circle()(vr_12q, tag );
    Bounded_side bs2 = Bounded_side_of_CCW_circle()(vr_1q2, tag );

    bool is_bs1 = (bs1 == ON_UNBOUNDED_SIDE);
    bool is_bs2 = (bs2 == ON_UNBOUNDED_SIDE);

    // both the ccw and cw circles do not exist
    if ( !is_bs1 && !is_bs2 ) {
      return b;
    }

    // the ccw circle exists but not the cw
    if ( is_bs1 && !is_bs2 ) {
      return b;
    }

    // the cw circle exists but not the ccw
    if ( !is_bs1 && is_bs2 ) {
      return b;
    }

    // both circles exist

    // check whether the shadow region is connected, i.e., wether it is
    // of the form (a, b) or (-oo, a) U (b, +oo)
    Bitangent_line bl_12(p1, p2);

    Sign stc =
      Sign_of_distance_from_bitangent_line()(bl_12, q, tag);

    Inverted_weighted_point u3 = inverter(p3);
    Bitangent_line blinv_23(u2, u3);

    CGAL_assertion( stc != ZERO );
    bool is_shadow_region_connected = (stc == POSITIVE);

    if ( is_shadow_region_connected ) {
      // the shadow region is of the form (a, b)
      if ( b ) { return false; }

      Voronoi_circle vc_123(blinv_23);
      Voronoi_circle vc_12q(vr_12q);

      Comparison_result r =
	Order_on_finite_bisector()(vc_123, vc_12q, p1, p2, tag);

      return ( r == SMALLER );
    }

    // the shadow region is of the form (-oo, a) U (b, +oo)
    if ( !b ) { return false; }

    Voronoi_circle vc_123(blinv_23);
    Voronoi_circle vc_1q2(vr_1q2);
    Comparison_result r =
      Order_on_finite_bisector()(vc_123, vc_1q2, p1, p2, tag);

    return ( r != SMALLER );
  }




  template<class Method_tag>
  bool
  operator()(const Weighted_point& p1,
	     const Weighted_point& p2,
	     const Weighted_point& q, bool b, Method_tag tag) {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::shadow_region_type_counter++;
#endif
    //
    Weighted_point_inverter inverter(p1);
    Inverted_weighted_point u2 = inverter(p2);
    Inverted_weighted_point v = inverter(q);

    Voronoi_radius vr_12q(u2, v);
    Voronoi_radius vr_1q2 = vr_12q.get_symmetric();

    Bounded_side bs1 = Bounded_side_of_CCW_circle()(vr_12q, tag );
    Bounded_side bs2 = Bounded_side_of_CCW_circle()(vr_1q2, tag );

    bool is_bs1 = (bs1 == ON_UNBOUNDED_SIDE);
    bool is_bs2 = (bs2 == ON_UNBOUNDED_SIDE);

    // both the ccw and cw circles do not exist
    if ( !is_bs1 && !is_bs2 ) {
      return b;
    }

    // only the ccw circle exists
    if ( is_bs1 && !is_bs2 ) { return false; }

    // only the cw circle exists
    if ( !is_bs1 && is_bs2 ) { return false; }

    // both circles exist
    
    // check whether the shadow region is connected, i.e., wether it is
    // of the form (a, b) or (-oo, a) U (b, +oo)

    return !b;
  }
};

//--------------------------------------------------------------------

template< class Point, class Weight >
class Bounded_side_of_CCW_circular_arc
{
public:
  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Weighted_point_inverter<Point, Weight>  Weighted_point_inverter;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef Voronoi_radius<Point, Weight>           Voronoi_radius;
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef typename Point::R::FT                          FT;
public:
  
  template< class Method_tag >
  Bounded_side operator()(const Bitangent_line& l1,
			  const Bitangent_line& l2,
			  const Bitangent_line& l3, Method_tag tag) const
    {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::inside_circular_arc_counter++;
#endif
      // This function checks whether the direction (a3, b3) (defined in
      // the unit circle) is inside the CCW circular arc defined on the
      // unit circle by the directions (a1, b1) and (a2, b2). By CCW arc
      // we mean that we walk on the unit circle in the CCW order from
      // (a1,b1) to (a2, b2) to define the arc.


      Orientation o = chi2(l1, l2, tag);//sign_of_determinant2x2(a1, b1, a2, b2);

      //      std::cout << "orientation(l1, l2): " << int(o) << std::endl;

      if ( o == COLLINEAR ) {
	Bitangent_line l2_rot = l2.get_rot90();
	Sign dot = chi2(l1, l2_rot, tag);//sign_of_determinant2x2(a1, b1, -b2, a2);
	CGAL_assertion( dot != ZERO );
	Orientation o1 = chi2(l1, l3, tag);//sign_of_determinant2x2(a1, b1, a3, b3);
	if ( dot == POSITIVE ) {
	  if ( o1 != COLLINEAR ) { return ON_UNBOUNDED_SIDE; }
	  Bitangent_line l3_rot = l3.get_rot90();
	  Sign dot1 = chi2(l1, l3_rot, tag);//sign_of_determinant2x2(a1, b1, -b3, a3);
	  CGAL_assertion( dot1 != ZERO );
	  return ( dot1 == POSITIVE ) ? ON_BOUNDARY : ON_UNBOUNDED_SIDE;
	}

	if ( o1 == LEFT_TURN ) { return ON_BOUNDED_SIDE; }
	return ( o1 == COLLINEAR ) ? ON_BOUNDARY : ON_UNBOUNDED_SIDE;
      } else if ( o == LEFT_TURN ) {
	Orientation o1 = chi2(l1, l3, tag);//sign_of_determinant2x2(a1, b1, a3, b3);
	Orientation o2 = chi2(l2, l3, tag);//sign_of_determinant2x2(a2, b2, a3, b3);

	//	std::cout << "orientation(l1, l3): " << int(o1) << std::endl;
	//	std::cout << "orientation(l2, l3): " << int(o2) << std::endl;

	if ( o1 == LEFT_TURN ) {
	  if ( o2 == COLLINEAR ) { return ON_BOUNDARY; }
	  return ( o2 == RIGHT_TURN ) ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
	} else if ( o1 == COLLINEAR ) {
	  CGAL_assertion( o2 != COLLINEAR );
	  return ( o2 == RIGHT_TURN ) ? ON_BOUNDARY : ON_UNBOUNDED_SIDE;
	}
	return ON_UNBOUNDED_SIDE;
      }
      Orientation o1 = chi2(l1, l3, tag);//sign_of_determinant2x2(a1, b1, a3, b3);
      Orientation o2 = chi2(l2, l3, tag);//sign_of_determinant2x2(a2, b2, a3, b3);

      //      std::cout << "orientation(l1, l3): " << int(o1) << std::endl;
      //      std::cout << "orientation(l2, l3): " << int(o2) << std::endl;


      if ( o1 == LEFT_TURN || o2 == RIGHT_TURN ) { return ON_BOUNDED_SIDE; }
      if ( o1 == COLLINEAR || o2 == COLLINEAR ) { return ON_BOUNDARY; }
      return ON_UNBOUNDED_SIDE;
    }

  Sign chi2(const Bitangent_line& bl1,
	    const Bitangent_line& bl2, Sqrt_field_tag) const
    {
      FT sigma = bl1.dx() * bl2.dx() + bl1.dy() * bl2.dy();
      FT delta = bl1.dx() * bl2.dy() - bl1.dy() * bl2.dx();

      //      FT E1 = -bl2.dw() * sigma;
      //      FT E2 = bl1.dw() * sigma;
      FT E1 = bl2.dw() * sigma;
      FT E2 = -bl1.dw() * sigma;
      FT E3 = delta;
      FT E4 = bl1.dw() * bl2.dw() * delta;
      FT p = bl1.delta();
      FT P = bl2.delta();

      FT E = E1 * CGAL_NTS sqrt(p) + E2 * CGAL_NTS sqrt(P)
	+ E3 * CGAL_NTS sqrt(p * P) + E4;

      return CGAL_NTS sign(E);
    }

  inline
  Sign chi2(const Bitangent_line& bl1,
	    const Bitangent_line& bl2, Ring_tag) const
    {
      return chi2(bl1.dx(), bl1.dy(), -bl1.dw(), bl1.d(), bl1.delta(),
		  bl2.dx(), bl2.dy(), -bl2.dw(), bl2.d(), bl2.delta());
    }

  Sign chi2(const FT& a, const FT& b, const FT& r,
	    const FT& d, const FT& p,
	    const FT& A, const FT& B, const FT& R,
	    const FT& D, const FT& P) const
    {
      FT sigma = a * A + b * B;
      FT delta = det2x2_by_formula(a, b, A, B);

      Sign sign_sigma = CGAL_NTS sign(sigma);
      Sign sign_delta = CGAL_NTS sign(delta);
      Sign sign_r = CGAL_NTS sign(r);
      Sign sign_R = CGAL_NTS sign(R);

      Sign sign_E1 = opposite(sign_R * sign_sigma);
      Sign sign_E2 = sign_r * sign_sigma;
      Sign sign_E3 = sign_delta;
      Sign sign_E4 = sign_r * sign_R * sign_delta;

      Sign sign_E1_plus_E3_P, sign_E4_plus_E2_P;

      //      FT d = CGAL_NTS square(a) + CGAL_NTS square(b);
      FT G = CGAL_NTS square(R) * d;
      FT delta2 = CGAL_NTS square(delta);

      if ( sign_E3 == ZERO ) {
	sign_E1_plus_E3_P = sign_E1;
      } else {
	if ( sign_E3 == sign_E1 ) {
	  sign_E1_plus_E3_P = sign_E3;
	} else {
	  FT F1 = delta2 - G;
	  sign_E1_plus_E3_P = sign_E3 * CGAL_NTS sign(F1);
	}
      }

      if ( sign_E2 == ZERO ) {
	sign_E4_plus_E2_P = sign_E4;
      } else {
	if ( sign_E2 == sign_E4 ) {
	  sign_E4_plus_E2_P = sign_E2;
	} else {
	  FT F2 = CGAL_NTS square(sigma) - G;
	  if ( sign_r == ZERO ) {
	    sign_E4_plus_E2_P = ZERO;
	  } else {
	    sign_E4_plus_E2_P = sign_E2 * CGAL_NTS sign(F2);
	  }
	}
      }

      if ( sign_E1_plus_E3_P == ZERO ) { return sign_E4_plus_E2_P; }
      if ( sign_E1_plus_E3_P == sign_E4_plus_E2_P ) {
	return sign_E1_plus_E3_P;
      }

      Sign sign_E5 = opposite(sign_R * sign_sigma * sign_delta);

      //      FT D = CGAL_NTS square(A) + CGAL_NTS square(B);
      //      FT P = D - CGAL_NTS square(R);

      FT F3 = P * delta2 + CGAL_NTS square(R * sigma) -
	CGAL_NTS square(r * D);

      Sign sign_E6 = CGAL_NTS sign(F3);

      if ( sign_E5 == ZERO ) { return sign_E1_plus_E3_P * sign_E6; }
      if ( sign_E5 == sign_E6 ) {
	return sign_E1_plus_E3_P * sign_E5;
      }

      //      FT p = d - CGAL_NTS square(r);
      FT rR = r * R;
      FT pP = p * P;
      //error();
      FT F4 = CGAL_NTS square(sigma - rR) - pP;
      FT F5 = CGAL_NTS square(sigma + rR) - pP;

      Sign sign_E7 = opposite(CGAL_NTS sign(F4) * CGAL_NTS sign(F5));

      return sign_E1_plus_E3_P * sign_E5 * sign_E7;
    }

};


//--------------------------------------------------------------------

template < class R, class W = typename R::RT >
class Infinite_edge_test
{
public:
  typedef Point_2<R>                              Point;
  typedef W                                       Weight;

  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Weighted_point_inverter<Point, Weight>  Weighted_point_inverter;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef Voronoi_radius<Point, Weight>           Voronoi_radius;
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef typename R::FT                          FT;

  typedef Bounded_side_of_CCW_circle<Point, Weight>
                                               Bounded_side_of_CCW_circle;
  typedef Sign_of_distance_from_bitangent_line<Point, Weight>
                                     Sign_of_distance_from_bitangent_line;
  typedef Sign_of_distance_from_CCW_circle<Point, Weight>
                                         Sign_of_distance_from_CCW_circle;
  typedef Order_on_finite_bisector<Point, Weight>
                                                 Order_on_finite_bisector;

  typedef Bounded_side_of_CCW_circular_arc<Point, Weight>
                                         Bounded_side_of_CCW_circular_arc;

public:
  template<class Method_tag>
  bool
  operator()(const Weighted_point& p2,
	     const Weighted_point& p3,
	     const Weighted_point& p4,
	     const Weighted_point& q, bool b, Method_tag tag) {
    Bitangent_line bl_32(p3, p2);
    Bitangent_line bl_24(p2, p4);
    Bitangent_line bl_2q(p2, q);

    Bounded_side bs1 =
      Bounded_side_of_CCW_circular_arc()(bl_32, bl_24, bl_2q, tag);

    if ( b ) {
      if ( bs1 == ON_BOUNDARY ) {
	Bitangent_line bl_q2(q, p2);
	Bounded_side bs2 =
	  Bounded_side_of_CCW_circular_arc()(bl_32, bl_24, bl_q2, tag);

	if ( bs2 != ON_BOUNDARY ) {
	  return ( bs2 != ON_BOUNDED_SIDE );
	}
	return !b;
      }

      return ( bs1 != ON_BOUNDED_SIDE );
    }
    if ( bs1 == ON_BOUNDARY ) {
      Bitangent_line bl_q2(q, p2);
      Bounded_side bs2 =
	Bounded_side_of_CCW_circular_arc()(bl_32, bl_24, bl_q2, tag);

      if ( bs2 != ON_BOUNDARY ) {
	return ( bs2 == ON_BOUNDED_SIDE );
      }
      return !b;
    }
    return ( bs1 == ON_BOUNDED_SIDE );
  }

};


//--------------------------------------------------------------------

template < class R, class W = typename R::RT >
class Is_degenerate_edge_test
{
public:
  typedef Point_2<R>                              Point;
  typedef W                                       Weight;


  typedef Weighted_point<Point, Weight>           Weighted_point;
  typedef Weighted_point_inverter<Point, Weight>  Weighted_point_inverter;
  typedef Inverted_weighted_point<Point, Weight>  Inverted_weighted_point;
  typedef Bitangent_line<Point, Weight>           Bitangent_line;
  typedef Voronoi_circle<Point, Weight>           Voronoi_circle;
  typedef typename R::FT                          FT;

  typedef Sign_of_distance_from_CCW_circle<Point, Weight>
                                          Sign_of_distance_from_CCW_circle;
  typedef Order_on_finite_bisector<Point, Weight>
                                                  Order_on_finite_bisector;

public:

  template<class Method_tag>
  inline bool operator()(const Weighted_point& p1,
			 const Weighted_point& p2,
			 const Weighted_point& p3,
			 const Weighted_point& p4,
			 Method_tag tag) const {
    Weighted_point_inverter inverter(p1);
    Inverted_weighted_point u2 = inverter(p2);
    Inverted_weighted_point u3 = inverter(p3);
    Inverted_weighted_point u4 = inverter(p4);

    Sign s;

    Bitangent_line blinv_23(u2, u3);
    s = Sign_of_distance_from_CCW_circle()(blinv_23, u4, tag);
    if ( s != ZERO ) { return false; }

    Bitangent_line blinv_42(u4, u2);
    s = Sign_of_distance_from_CCW_circle()(blinv_42, u3, tag);
    if ( s != ZERO ) { return false; }

    Voronoi_circle vc_123(blinv_23);
    Voronoi_circle vc_142(blinv_42);
    Comparison_result r =
      Order_on_finite_bisector()(vc_123, vc_142, p1, p2, tag);

    return ( r == EQUAL );
  }
};

//--------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_PREDICATES_C2_H

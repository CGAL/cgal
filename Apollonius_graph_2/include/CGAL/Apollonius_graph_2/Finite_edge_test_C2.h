// Copyright (c) 2003,2004,2006  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST_C2_H

#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Apollonius_graph_2/Predicate_constructions_C2.h>

#include <CGAL/functions_on_signs.h>
#include <CGAL/Apollonius_graph_2/compare_quadratic.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

//--------------------------------------------------------------------

template< class K >
class Orientation_wrt_symmetry_axis_2
{
public:
  typedef typename K::Point_2        Point_2;
  typedef Voronoi_circle_2<K>        Voronoi_circle;
  typedef typename K::FT             FT;
  typedef typename K::Orientation    Orientation;

public:

  Orientation
  operator()(const Voronoi_circle& vc, const Point_2& p1,
	     const Point_2& p2, const Field_with_sqrt_tag&) const
    {
      FT a = vc.a1() + vc.a2() * CGAL::sqrt(vc.delta());
      FT b = vc.b1() + vc.b2() * CGAL::sqrt(vc.delta());
      FT det = a * (p2.y() - p1.y()) - b * (p2.x() - p1.x());
      return CGAL::sign(det);
    }

  Orientation
  operator()(const Voronoi_circle& vc, const Point_2& p1,
	     const Point_2& p2, const Integral_domain_without_division_tag&) const
    {
      FT dx = p2.x() - p1.x();
      FT dy = p2.y() - p1.y();
      FT A = vc.a1() * dy - vc.b1() * dx;
      FT B = vc.a2() * dy - vc.b2() * dx;
      return sign_a_plus_b_x_sqrt_c(A, B, vc.delta());
    }
};


//--------------------------------------------------------------------

template< class K >
class Compare_Voronoi_radii_2
{
public:
  typedef Voronoi_circle_2<K>               Voronoi_circle;
  typedef typename K::FT                    FT;
  typedef typename K::Sign                  Sign;
  typedef typename K::Comparison_result     Comparison_result;

private:

  Sign sign_of_P4(const FT& u, const FT& v,
		  const FT& du, const FT& dv, const FT& dr,
		  const FT& Du, const FT& Dv, const FT& Dr) const
    {
      std::pair<FT,FT> factors = factors_of_P4(u, v, du, dv, dr,
					       Du, Dv, Dr);
      Sign s1 = CGAL::sign(factors.first);
      Sign s2 = CGAL::sign(factors.second);
      return s1 * s2;
    }

  std::pair<FT,FT>
  factors_of_P4(const FT& u, const FT& v,
		const FT& du, const FT& dv, const FT& dr,
		const FT& Du, const FT& Dv, const FT& Dr) const
    {
      FT u2 = CGAL::square(u);
      FT v2 = CGAL::square(v);

      FT du2 = CGAL::square(du);
      FT Du2 = CGAL::square(Du);

      FT dv2 = CGAL::square(dv);
      FT Dv2 = CGAL::square(Dv);

      FT dr2 = CGAL::square(dr);
      FT Dr2 = CGAL::square(Dr);

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
	- uU_P_vV * drDr - CGAL::square(uV_M_Uv);

      FT F2 = CGAL::square(u2_P_v2) * (du2_P_dv2 * Dr2 + Du2_P_Dv2 * dr2);
      F2 -= u2_P_v2 * (u2_M_v2 * uU_M_vV + uv * uV_P_Uv) * drDr;
      F2 -= CGAL::square(u2_M_v2 * uV_P_Uv - uv * uU_M_vV);

      std::pair<FT, FT> factors(F1,F2);
      return factors;
    }

public:
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     const Field_with_sqrt_tag&) const
    {
      FT c1 = (vc1.c1() + vc1.c2() * CGAL::sqrt(vc1.delta())) / vc1.d();
      FT c2 = (vc2.c1() + vc2.c2() * CGAL::sqrt(vc2.delta())) / vc2.d();

      Comparison_result r = CGAL::compare(c2, c1);
      return r;
    }

  // this is the naive way but without divisions and square roots; the
  // degree becomes 36 in this case.
  /*
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     Integral_domain_without_division_tag)
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
	     const Integral_domain_without_division_tag&) const
    {
      bool is_first_root1 = vc1.is_first_root();
      bool is_first_root2 = vc2.is_first_root();

      CGAL_precondition( CGAL::is_positive(vc1.alpha()) );
      CGAL_precondition( CGAL::is_positive(vc2.alpha()) );

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
	     Integral_domain_without_division_tag)
    {
      bool is_first_root1 = vc1.is_first_root();
      bool is_first_root2 = vc2.is_first_root();

      CGAL_precondition( CGAL::is_positive(vc1.alpha()) );
      CGAL_precondition( CGAL::is_positive(vc2.alpha()) );

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

template< class K >
class Order_on_finite_bisector_2
{
public:
  typedef Voronoi_circle_2<K>               Voronoi_circle;
  typedef typename K::Site_2                Site_2;
  typedef typename K::FT                    FT;
  typedef typename K::Comparison_result     Comparison_result;
  typedef typename K::Orientation           Orientation;

  typedef Compare_Voronoi_radii_2<K>        Compare_Voronoi_radii;

  typedef Orientation_wrt_symmetry_axis_2<K>
                                    Orientation_wrt_symmetry_axis;
  
public:
  template<class Method_tag>
  Comparison_result
  operator()(const Voronoi_circle& vc1, const Voronoi_circle& vc2,
	     const Site_2& p1, const Site_2& p2,
	     const Method_tag& tag) const
    {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::order_on_bisector_counter++;
#endif

      Orientation o1 =
	Orientation_wrt_symmetry_axis()(vc1, p1.point(), p2.point(), tag);
      Orientation o2 =
	Orientation_wrt_symmetry_axis()(vc2, p1.point(), p2.point(), tag);

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

template < class K >
class Finite_edge_interior_conflict
{
public:
  typedef typename K::Site_2                Site_2;
  typedef Weighted_point_inverter_2<K>      Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>      Inverted_weighted_point;
  typedef Voronoi_radius_2<K>               Voronoi_radius;
  typedef Voronoi_circle_2<K>               Voronoi_circle;
  typedef Bitangent_line_2<K>               Bitangent_line;
  typedef typename K::FT                    FT;
  typedef typename K::Sign                  Sign;
  typedef typename K::Bounded_side          Bounded_side;
  typedef typename K::Comparison_result     Comparison_result;

  typedef Bounded_side_of_CCW_circle_2<K>   Bounded_side_of_CCW_circle;
  typedef Order_on_finite_bisector_2<K>     Order_on_finite_bisector;

  typedef Sign_of_distance_from_bitangent_line_2<K>
                                     Sign_of_distance_from_bitangent_line;

  typedef Sign_of_distance_from_CCW_circle_2<K>
                                         Sign_of_distance_from_CCW_circle;

public:
  template<class Method_tag>
  bool
  operator()(const Site_2& p1,
	     const Site_2& p2,
	     const Site_2& p3,
	     const Site_2& p4,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
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

template < class K >
class Finite_edge_interior_conflict_degenerated
{
public:
  typedef typename K::Site_2                Site_2;
  typedef Weighted_point_inverter_2<K>      Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>      Inverted_weighted_point;
  typedef Voronoi_radius_2<K>               Voronoi_radius;
  typedef Voronoi_circle_2<K>               Voronoi_circle;
  typedef Bitangent_line_2<K>               Bitangent_line;
  typedef typename K::FT                    FT;
  typedef typename K::Sign                  Sign;
  typedef typename K::Comparison_result     Comparison_result;
  typedef typename K::Bounded_side          Bounded_side;

  typedef Bounded_side_of_CCW_circle_2<K>   Bounded_side_of_CCW_circle;
  typedef Order_on_finite_bisector_2<K>     Order_on_finite_bisector;

  typedef Sign_of_distance_from_bitangent_line_2<K>
                                     Sign_of_distance_from_bitangent_line;

  typedef Sign_of_distance_from_CCW_circle_2<K>
                                         Sign_of_distance_from_CCW_circle;
public:
  template<class Method_tag>
  bool
  operator()(const Site_2& p1, const Site_2& p2, const Site_2& p3,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
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

    Sign stc = Sign_of_distance_from_bitangent_line()(bl_12, q, tag);

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
  operator()(const Site_2& p1, const Site_2& p2,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
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


template<class K, class MTag>
class Finite_edge_interior_conflict_2
{
public:
  typedef K                      Kernel;
  typedef MTag                   Method_tag;

  typedef typename K::Site_2     Site_2;

private:
  typedef Finite_edge_interior_conflict_degenerated<K>   Test_degenerated;
  typedef Finite_edge_interior_conflict<K>               Test;

public:
  typedef bool                  result_type;
  struct argument_type {};

  inline
  bool operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& q, bool b) const
  {
    return Test_degenerated()(p1, p2, q, b, Method_tag());
  }

  inline
  bool operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& p3, const Site_2& q, bool b) const
  {
    return Test_degenerated()(p1, p2, p3, q, b, Method_tag());    
  }

  inline
  bool operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& p3, const Site_2& p4,
		  const Site_2& q, bool b) const
  {
    return Test()(p1, p2, p3, p4, q, b, Method_tag());    
  }
};


//--------------------------------------------------------------------

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST_C2_H

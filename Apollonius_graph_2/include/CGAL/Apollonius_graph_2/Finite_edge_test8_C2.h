// Copyright (c) 2007 INRIA Sophia-Antipolis (France).
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


#ifndef CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST8_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST8_C2_H

#include <CGAL/Apollonius_graph_2/Finite_edge_test_C2.h>
#include <CGAL/Apollonius_graph_2/Orientation8_C2.h>
#include <CGAL/Apollonius_graph_2/Oriented_side_of_bisector_C2.h>


namespace CGAL {

namespace ApolloniusGraph_2 {

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template < class K, class MTag >
class Inside_Voronoi_quadrilateral8
{
public:
  typedef K                                 Kernel;
  typedef MTag                              Method_tag;
  typedef typename K::Site_2                Site_2;
  typedef typename K::Orientation           Orientation;
  typedef typename K::Oriented_side         Oriented_side;

private:
  typedef Orientation8_C2<K,MTag>           Orientation_2;
  typedef Constructive_orientation8_C2<K,MTag>   Constructive_orientation_2;
  typedef Oriented_side_of_bisector_2<K,MTag>    Side_of_bisector_2;

public:
  bool
  operator()(const Site_2& p1,
	     const Site_2& p2,
	     const Site_2& p3,
	     const Site_2& p4,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
#if 1
    Constructive_orientation_2 orientation123(p1, p2, p3, true);
    Constructive_orientation_2 orientation142(p1, p4, p2, false);

    Orientation o123_s = orientation123();
    Orientation o142_s = orientation142();

#if 0
#ifndef NDEBUG
    Orientation o123_s1 = orientation123(p1, p2);
    Orientation o142_s1 = orientation142(p1, p2);

    CGAL_assertion( o123_s == o123_s1 );
    CGAL_assertion( o142_s == o142_s1 );
#endif
#endif

    // first we consider the case where both Voronoi circles are in
    // conflict; we want to determine if the entire edge is in
    // conflict

    if ( b ) {
      if ( o123_s != o142_s ) {
	Orientation o123_1 = orientation123(p1, q);
	Orientation o123_2 = orientation123(p2, q);
	if ( o123_1 == POSITIVE && o123_2 == NEGATIVE ) { return true; }

	Orientation o142_1 = orientation142(p1, q);
	Orientation o142_2 = orientation142(p2, q);
	return (o142_1 == NEGATIVE && o142_2 == POSITIVE);
      }

      Oriented_side os = Side_of_bisector_2()(p1, p2, q.point());

      if ( o123_s == POSITIVE ) {
	if ( os == ON_POSITIVE_SIDE ) {
	  Orientation o142_1 = orientation142(p1, q);
	  if ( o142_1 == NEGATIVE ) { return true; }

	  Orientation o123_1 = orientation123(p1, q);
	  return o123_1 == POSITIVE;
	}

	Orientation o142_2 = orientation142(p2, q);
	if ( o142_2 == POSITIVE ) { return true; }

	Orientation o123_2 = orientation123(p2, q);
	return o123_2 == NEGATIVE;
      }

      if ( o123_s == NEGATIVE ) {
	if ( os == ON_POSITIVE_SIDE ) {
	  Orientation o123_1 = orientation123(p1, q);
	  if ( o123_1 == POSITIVE ) { return true; }

	  Orientation o142_1 = orientation142(p1, q);
	  return o142_1 == NEGATIVE;
	}

	Orientation o123_2 = orientation123(p2, q);	
	if ( o123_2 == NEGATIVE ) { return true; }

	Orientation o142_2 = orientation142(p2, q);
	return o142_2 == POSITIVE;
      }

      CGAL_assertion( o123_s == ZERO );
      return true;
    }

    // now consider the case where the two Voronoi circles are not in
    // conflict; we want to determine if part of the interior is in
    // conflict

    if ( o123_s != o142_s ) {
      Orientation o123_1 = orientation123(p1, q);
      if ( o123_1 != POSITIVE ) { return false; }

      Orientation o123_2 = orientation123(p2, q);
      if ( o123_2 != NEGATIVE ) { return false; }
      
      Orientation o142_1 = orientation142(p1, q);
      if ( o142_1 != NEGATIVE ) { return false; }

      Orientation o142_2 = orientation142(p2, q);
      return o142_2 == POSITIVE;
    }

    Oriented_side os = Side_of_bisector_2()(p1, p2, q.point());

    if ( o123_s == POSITIVE ) {
      if ( os == ON_POSITIVE_SIDE ) {
	Orientation o123_1 = orientation123(p1, q);
	if ( o123_1 != POSITIVE ) { return false; }

	Orientation o142_1 = orientation142(p1, q);
	return o142_1 == NEGATIVE;
      }

      Orientation o123_2 = orientation123(p2, q);
      if ( o123_2 != NEGATIVE ) { return false; }

      Orientation o142_2 = orientation142(p2, q);
      return o142_2 == POSITIVE;
    }

    if ( o123_s == NEGATIVE ) {
      if ( os == ON_POSITIVE_SIDE ) {
	Orientation o142_1 = orientation142(p1, q);
	if ( o142_1 != NEGATIVE ) { return false; }

	Orientation o123_1 = orientation123(p1, q);
	return o123_1 == POSITIVE;
      }
      Orientation o142_2 = orientation142(p2, q);
      if ( o142_2 != POSITIVE ) { return false; }

      Orientation o123_2 = orientation123(p2, q);
      return o123_2 == NEGATIVE;
    }

    CGAL_assertion( o123_s == ZERO );
    return false;

#if 0
    if ( os == ON_POSITIVE_SIDE ) {
      Orientation o123_1 = orientation(p1, p2, p3, p1, q);
      Orientation o142_1 = orientation(p1, p4, p2, p1, q);
      //      std::cerr << "o123_1: " << o123_1 << std::endl;
      //      std::cerr << "o142_1: " << o142_1 << std::endl;

      if ( b ) {
	return !(o123_1 == NEGATIVE && o142_1 == POSITIVE);
      }
      return o123_1 == POSITIVE && o142_1 == NEGATIVE;
    }

    Orientation o123_2 = orientation(p1, p2, p3, p2, q);
    Orientation o142_2 = orientation(p1, p4, p2, p2, q);
    //    std::cerr << "o123_2: " << o123_2 << std::endl;
    //    std::cerr << "o142_2: " << o142_2 << std::endl;

    if ( b ) {
      return !(o123_2 == POSITIVE && o142_2 == NEGATIVE);
    }
    return o123_2 == NEGATIVE && o142_2 == POSITIVE;
#endif
#else
    Constructive_orientation_2 orientation123(p1, p2, p3);
    Constructive_orientation_2 orientation142(p1, p4, p2);

    Orientation o123_s = orientation123(p1, p2);
    Orientation o142_s = orientation142(p1, p2);

    // first we consider the case where both Voronoi circles are in
    // conflict; we want to determine if the entire edge is in
    // conflict

    if ( b ) {
      if ( o123_s != o142_s ) {
	Orientation o123_1 = orientation123(p1, q);
	Orientation o123_2 = orientation123(p2, q);
	if ( o123_1 == POSITIVE && o123_2 == NEGATIVE ) { return true; }

	Orientation o142_1 = orientation142(p1, q);
	Orientation o142_2 = orientation142(p2, q);
	return (o142_1 == NEGATIVE && o142_2 == POSITIVE);
      }

      if ( o123_s == POSITIVE ) {
	Orientation o142_1 = orientation142(p1, q);
	if ( o142_1 == NEGATIVE ) { return true; }

	Orientation o142_2 = orientation142(p2, q);
	if ( o142_2 == POSITIVE ) { return true; }

	Orientation o123_1 = orientation123(p1, q);
	if ( o123_1 != POSITIVE ) { return false; }

	Orientation o123_2 = orientation123(p2, q);
	return o123_2 == NEGATIVE;
      }

      if ( o123_s == NEGATIVE ) {
	Orientation o123_1 = orientation123(p1, q);
	if ( o123_1 == POSITIVE ) { return true; }

	Orientation o123_2 = orientation123(p2, q);
	if ( o123_2 == NEGATIVE ) { return true; }

	Orientation o142_1 = orientation142(p1, q);
	if ( o142_1 != NEGATIVE ) { return false; }

	Orientation o142_2 = orientation142(p2, q);
	return o142_2 == POSITIVE;
      }

      CGAL_assertion( o123_s == ZERO );
      return true;
    }

    // now consider the case where the two Voronoi circles are not in
    // conflict; we want to determine if part of the interior is in
    // conflict

    if ( o123_s != o142_s ) {
      Orientation o123_1 = orientation123(p1, q);
      if ( o123_1 != POSITIVE ) { return false; }

      Orientation o123_2 = orientation123(p2, q);
      if ( o123_2 != NEGATIVE ) { return false; }
      
      Orientation o142_1 = orientation142(p1, q);
      if ( o142_1 != NEGATIVE ) { return false; }

      Orientation o142_2 = orientation142(p2, q);
      return o142_2 == POSITIVE;
    }

    if ( o123_s == POSITIVE ) {
      Orientation o123_1 = orientation123(p1, q);
      if ( o123_1 != POSITIVE ) { return false; }

      Orientation o123_2 = orientation123(p2, q);
      if ( o123_2 != NEGATIVE ) { return false; }

      Orientation o142_1 = orientation142(p1, q);
      if ( o142_1 == NEGATIVE ) { return true; }

      Orientation o142_2 = orientation142(p2, q);
      return o142_2 == POSITIVE;
    }

    if ( o123_s == NEGATIVE ) {
      Orientation o142_1 = orientation142(p1, q);
      if ( o142_1 != NEGATIVE ) { return false; }

      Orientation o142_2 = orientation142(p2, q);
      if ( o142_2 != POSITIVE ) { return false; }

      Orientation o123_1 = orientation123(p1, q);
      if ( o123_1 == POSITIVE ) { return true; }

      Orientation o123_2 = orientation123(p2, q);
      return o123_2 == NEGATIVE;
    }

    CGAL_assertion( o123_s == ZERO );
    return false;
#endif
  }
};

//--------------------------------------------------------------------

template < class K, class MTag >
class Finite_edge_interior_conflict8
{
public:
  typedef MTag                              Method_tag;

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

  typedef Sign_of_distance_from_bitangent_line_2<K>
                                     Sign_of_distance_from_bitangent_line;

  typedef Sign_of_distance_from_CCW_circle_2<K>
                                         Sign_of_distance_from_CCW_circle;

public:
  bool
  operator()(const Site_2& p1,
	     const Site_2& p2,
	     const Site_2& p3,
	     const Site_2& p4,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
    typedef Inside_Voronoi_quadrilateral8<K,Method_tag> Inside_quadrilateral;

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

#if 0
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
#else
      return Inside_quadrilateral()(p1, p2, p3, p4, q, b, tag);
#endif
    }

    // the shadow region is of the form (-oo, a) U (b, +oo)
    if ( !b ) { return false; }

#if 0
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
#else
    return Inside_quadrilateral()(p1, p2, p3, p4, q, b, tag);
#endif
  }
};

//--------------------------------------------------------------------

template < class K, class MTag >
class Finite_edge_interior_conflict_degenerated8
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
  typedef typename K::Orientation           Orientation;
  typedef typename K::Bounded_side          Bounded_side;
  typedef typename K::Oriented_side         Oriented_side;

  typedef Bounded_side_of_CCW_circle_2<K>   Bounded_side_of_CCW_circle;
  typedef Order_on_finite_bisector_2<K>     Order_on_finite_bisector;
  typedef Orientation8_C2<K,MTag>           Orientation_2;

  typedef Sign_of_distance_from_bitangent_line_2<K>
                                     Sign_of_distance_from_bitangent_line;

  typedef Sign_of_distance_from_CCW_circle_2<K>
                                         Sign_of_distance_from_CCW_circle;

  typedef Oriented_side_of_bisector_2<K,MTag>    Side_of_bisector_2;

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

  inline
  Sign sqrt_ext_sign(const FT& A1, const FT& Dr, const FT& B,
		     const FT& Dx1, const FT& Dy1, 
		     const FT& Dx2, const FT& Dy2, 
		     const Field_with_sqrt_tag&) const
  {
    FT D = CGAL::square(Dx1) + CGAL::square(Dy1) - CGAL::square(Dr);
    return CGAL::sign(A1 * Dr + B * CGAL::sqrt(D));
  }

  inline
  Sign sqrt_ext_sign(const FT& A1, const FT& Dr, const FT& B,
		     const FT& Dx1, const FT& Dy1, 
		     const FT& Dx2, const FT& Dy2,
		     const Integral_domain_without_division_tag&) const
  {
    Sign sA = CGAL::sign(A1) * CGAL::sign(Dr);
    Sign sB = CGAL::sign(B);

    if ( sA == CGAL::ZERO ) { return sB; }
    if ( sB == CGAL::ZERO ) { return sA; }
    if ( sA == sB ) { return sA; }

    FT C = (CGAL::square(Dx2) + CGAL::square(Dy2)) * CGAL::square(Dr);
    return sA * CGAL::sign(C - CGAL::square(B));
  }

  template<class Method_tag>
  Orientation
  orientation_wrt_bitangent_perp(const Site_2& p1, const Site_2& p2,
				 const Site_2& q1, const Site_2& q2,
				 const Method_tag& tag) const
  {
    // computes the orientation predicate of the line perpendicular to
    // the bitangent of p1 and p2, passing through the center of q1,
    // and the center of q2

    FT Dx1 = p1.x() - p2.x();
    FT Dy1 = p1.y() - p2.y();
    FT Dr = p1.weight() - p2.weight();

    FT Dx2 = q1.x() - q2.x();
    FT Dy2 = q1.y() - q2.y();

    FT A1 = Dx1 * Dy2 - Dy1 * Dx2;
    FT B = Dx1 * Dx2 + Dy1 * Dy2;

    return sqrt_ext_sign(A1, Dr, B, Dx1, Dy1, Dx2, Dy2, tag);
  }



  template<class Method_tag>
  bool
  operator()(const Site_2& p1, const Site_2& p2,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
#if 0
    Side_of_bisector_2 side_of_bisector;
    Orientation o12_sym = Orientation_2()(p1, p2, q);

    Oriented_side os = side_of_bisector(p1, p2, q.point());

    if ( o12_sym == CGAL::POSITIVE ) {
      if ( os == ON_POSITIVE_SIDE ) {
	Orientation o21_1 = orientation_wrt_bitangent_perp(p2, p1, p1, q, tag);
	return o21_1 == NEGATIVE;
      }

      Orientation o21_2 = orientation_wrt_bitangent_perp(p2, p1, p2, q, tag);
      return o21_2 == POSITIVE;

    } else {
      if ( os == ON_POSITIVE_SIDE ) {
	Orientation o12_1 = orientation_wrt_bitangent_perp(p1, p2, p1, q, tag);
	return o12_1 == POSITIVE;
      }

      Orientation o12_2 = orientation_wrt_bitangent_perp(p1, p2, p2, q, tag);
      return o12_2 == NEGATIVE;
    }
#else
    Orientation o12_sym = Orientation_2()(p1, p2, q);

    if ( o12_sym == CGAL::POSITIVE ) {
      Orientation o21_1 = orientation_wrt_bitangent_perp(p2, p1, p1, q, tag);

      if ( o21_1 != NEGATIVE ) { return false; }

      Orientation o21_2 = orientation_wrt_bitangent_perp(p2, p1, p2, q, tag);

      return o21_2 == POSITIVE;
    } else {
      Orientation o12_1 = orientation_wrt_bitangent_perp(p1, p2, p1, q, tag);

      if ( o12_1 != POSITIVE ) { return false; }

      Orientation o12_2 = orientation_wrt_bitangent_perp(p1, p2, p2, q, tag);

      return o12_2 == NEGATIVE;
    }
#endif
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::shadow_region_type_counter++;
#endif
  }
};

//--------------------------------------------------------------------


template<class K, class MTag>
class Finite_edge_interior_conflict8_2
{
public:
  typedef K                      Kernel;
  typedef MTag                   Method_tag;

  typedef typename K::Site_2     Site_2;

private:
  typedef Finite_edge_interior_conflict_degenerated8<K,MTag>  Test_degenerated;
  typedef Finite_edge_interior_conflict8<K,MTag>              Test;

public:
  typedef bool                  result_type;
  struct argument_type {};

  inline
  bool operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& q, bool b) const
  {
    typedef Finite_edge_interior_conflict_2<K,MTag> Old_Test;
    bool t = Test_degenerated()(p1, p2, q, b, Method_tag());    
    //    bool t = Old_Test()(p1, p2, q, b);
#ifndef NDEBUG
    bool t_old = Old_Test()(p1, p2, q, b);

    if ( t != t_old ) {
      std::cerr << std::endl;
      std::cerr << "b: " << b << "; t: " << t
		<< "; t_old: " << t_old << std::endl;
      std::cerr << "p1: " << p1 << std::endl;
      std::cerr << "p2: " << p2 << std::endl;
      std::cerr << "q: " << q << std::endl;
    }

    CGAL_assertion( t == t_old );
#endif
    return t;
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
    typedef Finite_edge_interior_conflict_2<K,MTag> Old_Test;
    bool t = Test()(p1, p2, p3, p4, q, b, Method_tag());    
    //    bool t = Old_Test()(p1, p2, p3, p4, q, b);
#ifndef NDEBUG
    bool t_old = Old_Test()(p1, p2, p3, p4, q, b);

    if ( t != t_old ) {
      Oriented_side_of_bisector_2<K,MTag> side_of_bisector;

      Oriented_side os = side_of_bisector(p1, p2, q.point());

      std::cerr << "b: " << b << std::endl;
      std::cerr << "t: " << t << "; t_old: " << t_old << std::endl;
      std::cerr << "os: " << os << std::endl;
      std::cerr << "p1: " << p1 << std::endl;
      std::cerr << "p2: " << p2 << std::endl;
      std::cerr << "p3: " << p3 << std::endl;
      std::cerr << "p4: " << p4 << std::endl;
      std::cerr << "q: " << q << std::endl;
    }

    CGAL_assertion( t == t_old );
#endif
    return t; //Test()(p1, p2, p3, p4, q, b, Method_tag());    
  }
};


//--------------------------------------------------------------------

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST8_C2_H

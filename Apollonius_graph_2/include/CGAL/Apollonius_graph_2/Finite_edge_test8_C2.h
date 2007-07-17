// Copyright (c) 2007 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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


CGAL_BEGIN_NAMESPACE

CGAL_APOLLONIUS_GRAPH_2_BEGIN_NAMESPACE

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template < class K, class MTag >
class Finite_edge_interior_conflict8
{
public:
  typedef K                                 Kernel;
  typedef MTag                              Method_tag;
  typedef typename K::Site_2                Site_2;
  //  typedef typename K::FT                    FT;
  //  typedef typename K::Sign                  Sign;
  typedef typename K::Orientation           Orientation;
  //  typedef typename K::Comparison_result     Comparison_result;

  typedef Orientation8_C2<K,MTag>           Orientation_2;

public:
  bool
  operator()(const Site_2& p1,
	     const Site_2& p2,
	     const Site_2& p3,
	     const Site_2& p4,
	     const Site_2& q, bool b, const Method_tag& tag) const
  {
    Orientation_2 orientation;

    Orientation o123_s = orientation(p1, p2, p3, p1, p2);
    Orientation o142_s = orientation(p1, p4, p2, p1, p2);

    Orientation o123_1 = orientation(p1, p2, p3, p1, q);
    Orientation o142_1 = orientation(p1, p4, p2, p1, q);
    Orientation o123_2 = orientation(p1, p2, p3, p2, q);
    Orientation o142_2 = orientation(p1, p4, p2, p2, q);

    // first we consider the case where both Voronoi circles are in
    // conflict; we want to determine if the entire edge is in
    // conflict

    if ( b ) {
      if ( o123_s != o142_s ) {
	return !((o123_1 == NEGATIVE || o123_2 == POSITIVE) &&
		 (o142_1 == POSITIVE || o142_2 == NEGATIVE));
      }

      if ( o123_s == POSITIVE ) {
	return !((o123_1 == NEGATIVE || o123_2 == POSITIVE) &&
		 (o142_1 == POSITIVE && o142_2 == NEGATIVE));
      }

      if ( o123_s == NEGATIVE ) {
	return !((o123_1 == NEGATIVE && o123_2 == POSITIVE) &&
		 (o142_1 == POSITIVE || o142_2 == NEGATIVE));
      }

      CGAL_assertion( o123_s == ZERO );
      return true;      
    }

    // now consider the case where the two Voronoi circles are not in
    // conflict; we want to determine if part of the interior is in
    // conflict

    if ( o123_s != o142_s ) {
      return ((o123_1 == POSITIVE && o123_2 == NEGATIVE) &&
	      (o142_1 == NEGATIVE && o142_2 == POSITIVE));
    }

    if ( o123_s == POSITIVE ) {
      return ((o123_1 == POSITIVE && o123_2 == NEGATIVE) &&
	      (o142_1 == NEGATIVE || o142_2 == POSITIVE));
    }

    if ( o123_s == NEGATIVE ) {
      return ((o123_1 == POSITIVE || o123_2 == NEGATIVE) &&
	      (o142_1 == NEGATIVE && o142_2 == POSITIVE));
    }

    CGAL_assertion( o123_s == ZERO );
    return false;
  }
};

//--------------------------------------------------------------------

template < class K >
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
class Finite_edge_interior_conflict8_2
{
public:
  typedef K                      Kernel;
  typedef MTag                   Method_tag;

  typedef typename K::Site_2     Site_2;

private:
  typedef Finite_edge_interior_conflict_degenerated8<K>   Test_degenerated;
  typedef Finite_edge_interior_conflict8<K,MTag>          Test;

public:
  typedef bool                  result_type;
  struct argument_type {};
  struct Arity {};

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
    typedef Finite_edge_interior_conflict_2<K,MTag> Old_Test;
    bool t = Test()(p1, p2, p3, p4, q, b, Method_tag());    
    bool t_old = Old_Test()(p1, p2, p3, p4, q, b);

    if ( t != t_old ) {
      std::cerr << "t: " << t << "; t_old: " << t_old << std::endl;
      std::cerr << "p1: " << p1 << std::endl;
      std::cerr << "p2: " << p2 << std::endl;
      std::cerr << "p3: " << p3 << std::endl;
      std::cerr << "p4: " << p4 << std::endl;
      std::cerr << "q: " << q << std::endl;
    }

    CGAL_assertion( t == t_old );

    return t; //Test()(p1, p2, p3, p4, q, b, Method_tag());    
  }
};


//--------------------------------------------------------------------

CGAL_APOLLONIUS_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_TEST8_C2_H

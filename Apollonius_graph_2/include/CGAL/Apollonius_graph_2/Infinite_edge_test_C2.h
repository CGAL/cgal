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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_TEST_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_TEST_C2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Apollonius_graph_2/Predicate_constructions_C2.h>


namespace CGAL {

namespace ApolloniusGraph_2 {

//--------------------------------------------------------------------

template< class K >
class Bounded_side_of_CCW_circular_arc_2
{
public:
  typedef Weighted_point_inverter_2<K>    Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>    Inverted_weighted_point;
  typedef Voronoi_radius_2<K>             Voronoi_radius;
  typedef Voronoi_circle_2<K>             Voronoi_circle;
  typedef Bitangent_line_2<K>             Bitangent_line;
  typedef typename K::FT                  FT;
  typedef typename K::Bounded_side        Bounded_side;
  typedef typename K::Orientation         Orientation;
  typedef typename K::Sign                Sign;

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


      Orientation o = chi2(l1, l2, tag);
                   //sign_of_determinant(a1, b1, a2, b2);


      if ( o == COLLINEAR ) {
	Bitangent_line l2_rot = l2.get_rot90();
	Sign dot = chi2(l1, l2_rot, tag);
                              //sign_of_determinant(a1, b1, -b2, a2);

	CGAL_assertion( dot != ZERO );
	Orientation o1 = chi2(l1, l3, tag);
                          //sign_of_determinant(a1, b1, a3, b3);

	if ( dot == POSITIVE ) {
	  if ( o1 != COLLINEAR ) { return ON_UNBOUNDED_SIDE; }
	  Bitangent_line l3_rot = l3.get_rot90();
	  Sign dot1 = chi2(l1, l3_rot, tag);
                       //sign_of_determinant(a1, b1, -b3, a3);

	  CGAL_assertion( dot1 != ZERO );
	  return ( dot1 == POSITIVE ) ? ON_BOUNDARY : ON_UNBOUNDED_SIDE;
	}

	if ( o1 == LEFT_TURN ) { return ON_BOUNDED_SIDE; }
	return ( o1 == COLLINEAR ) ? ON_BOUNDARY : ON_UNBOUNDED_SIDE;
      } else if ( o == LEFT_TURN ) {
	Orientation o1 = chi2(l1, l3, tag);
	                //sign_of_determinant(a1, b1, a3, b3);
	Orientation o2 = chi2(l2, l3, tag);
                        //sign_of_determinant(a2, b2, a3, b3);

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
      Orientation o1 = chi2(l1, l3, tag);
                        //sign_of_determinant(a1, b1, a3, b3);
      Orientation o2 = chi2(l2, l3, tag);
                        //sign_of_determinant(a2, b2, a3, b3);

      //      std::cout << "orientation(l1, l3): " << int(o1) << std::endl;
      //      std::cout << "orientation(l2, l3): " << int(o2) << std::endl;


      if ( o1 == LEFT_TURN || o2 == RIGHT_TURN ) { return ON_BOUNDED_SIDE; }
      if ( o1 == COLLINEAR || o2 == COLLINEAR ) { return ON_BOUNDARY; }
      return ON_UNBOUNDED_SIDE;
    }

  Sign chi2(const Bitangent_line& bl1,
	    const Bitangent_line& bl2, Field_with_sqrt_tag) const
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

      FT E = E1 * CGAL::sqrt(p) + E2 * CGAL::sqrt(P)
	+ E3 * CGAL::sqrt(p * P) + E4;

      return CGAL::sign(E);
    }

  inline
  Sign chi2(const Bitangent_line& bl1,
	    const Bitangent_line& bl2, Integral_domain_without_division_tag) const
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
      FT delta = determinant(a, b, A, B);

      Sign sign_sigma = CGAL::sign(sigma);
      Sign sign_delta = CGAL::sign(delta);
      Sign sign_r = CGAL::sign(r);
      Sign sign_R = CGAL::sign(R);

      Sign sign_E1 = -sign_R * sign_sigma;
      Sign sign_E2 = sign_r * sign_sigma;
      Sign sign_E3 = sign_delta;
      Sign sign_E4 = sign_r * sign_R * sign_delta;

      Sign sign_E1_plus_E3_P, sign_E4_plus_E2_P;

      //      FT d = CGAL::square(a) + CGAL::square(b);
      FT G = CGAL::square(R) * d;
      FT delta2 = CGAL::square(delta);

      if ( sign_E3 == ZERO ) {
	sign_E1_plus_E3_P = sign_E1;
      } else {
	if ( sign_E3 == sign_E1 ) {
	  sign_E1_plus_E3_P = sign_E3;
	} else {
	  FT F1 = delta2 - G;
	  sign_E1_plus_E3_P = sign_E3 * CGAL::sign(F1);
	}
      }

      if ( sign_E2 == ZERO ) {
	sign_E4_plus_E2_P = sign_E4;
      } else {
	if ( sign_E2 == sign_E4 ) {
	  sign_E4_plus_E2_P = sign_E2;
	} else {
	  FT F2 = CGAL::square(sigma) - G;
	  if ( sign_r == ZERO ) {
	    sign_E4_plus_E2_P = ZERO;
	  } else {
	    sign_E4_plus_E2_P = sign_E2 * CGAL::sign(F2);
	  }
	}
      }

      if ( sign_E1_plus_E3_P == ZERO ) { return sign_E4_plus_E2_P; }
      if ( sign_E1_plus_E3_P == sign_E4_plus_E2_P ) {
	return sign_E1_plus_E3_P;
      }

      Sign sign_E5 = -sign_R * sign_sigma * sign_delta;

      //      FT D = CGAL::square(A) + CGAL::square(B);
      //      FT P = D - CGAL::square(R);

      FT F3 = P * delta2 + CGAL::square(R * sigma) - CGAL::square(r * D);

      Sign sign_E6 = CGAL::sign(F3);

      if ( sign_E5 == ZERO ) {
	return sign_E1_plus_E3_P * sign_E6;
      }
      if ( sign_E5 == sign_E6 ) {
	return sign_E1_plus_E3_P * sign_E5;
      }

      //      FT p = d - CGAL::square(r);
      FT rR = r * R;
      FT pP = p * P;
      //error();
      FT F4 = CGAL::square(sigma - rR) - pP;
      FT F5 = CGAL::square(sigma + rR) - pP;

      Sign sign_E7 = -CGAL::sign(F4) * CGAL::sign(F5);

      return sign_E1_plus_E3_P * sign_E5 * sign_E7;
    }

};


//--------------------------------------------------------------------

template < class K, class MTag >
class Infinite_edge_interior_conflict_2
{
public:
  typedef K                                Kernel;
  typedef MTag                             Method_tag;

  typedef typename K::Site_2               Site_2;
  typedef Weighted_point_inverter_2<K>     Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>     Inverted_weighted_point;
  typedef Voronoi_radius_2<K>              Voronoi_radius;
  typedef Voronoi_circle_2<K>              Voronoi_circle;
  typedef Bitangent_line_2<K>              Bitangent_line;
  typedef typename K::FT                   FT;
  typedef typename K::Bounded_side         Bounded_side;

  //  typedef CGAL::Bounded_side_of_CCW_circle<K>  Bounded_side_of_CCW_circle;
  //  typedef CGAL::Sign_of_distance_from_bitangent_line<K>
  //                                     Sign_of_distance_from_bitangent_line;
  //  typedef CGAL::Sign_of_distance_from_CCW_circle<K>
  //                                         Sign_of_distance_from_CCW_circle;
  //  typedef CGAL::Order_on_finite_bisector<K>      Order_on_finite_bisector;

  typedef Bounded_side_of_CCW_circular_arc_2<K>
                                         Bounded_side_of_CCW_circular_arc;

public:
  typedef bool            result_type;
  struct argument_type {};


  bool operator()(const Site_2& p2, const Site_2& p3,
		  const Site_2& p4, const Site_2& q, bool b) const
  {
    Method_tag tag = Method_tag();
    
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

} //namespace ApolloniusGraph_2

} //namespace CGAL


#endif // CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_TEST_C2_H

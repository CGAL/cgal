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



#ifndef CGAL_APOLLONIUS_GRAPH_2_INCIRCLE_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_INCIRCLE_C2_H

#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Apollonius_graph_2/Predicate_constructions_C2.h>

#include <CGAL/Apollonius_graph_2/Bounded_side_of_ccw_circle_C2.h>

#include <CGAL/functions_on_signs.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

//--------------------------------------------------------------------

template< class K >
class Sign_of_distance_from_bitangent_line_2
{
public:
  typedef Bitangent_line_2<K>               Bitangent_line;
  typedef typename K::Site_2                Site_2;
  typedef Inverted_weighted_point_2<K>      Inverted_weighted_point;
  typedef typename K::FT                    FT;
  typedef typename K::Sign                  Sign;

public:

  inline Sign
  operator()(const Bitangent_line& bl, const Site_2& q,
	     const Field_with_sqrt_tag&) const
    {
#ifdef AG2_PROFILE_PREDICATES
      ag2_predicate_profiler::distance_from_bitangent_counter++;
#endif
      FT a = bl.a1() + bl.a2() * CGAL::sqrt(bl.delta());
      FT b = bl.b1() + bl.b2() * CGAL::sqrt(bl.delta());
      FT c = bl.c1() + bl.c2() * CGAL::sqrt(bl.delta());
      FT r = a * q.x() + b * q.y() + c - q.weight() * bl.d();
      return CGAL::sign(r);
    }

  inline Sign
  operator()(const Bitangent_line& bl, const Site_2& q,
	     const Integral_domain_without_division_tag&) const
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


template< class K >
class Sign_of_distance_from_CCW_circle_2
{
public:
  typedef Bitangent_line_2<K>            Bitangent_line;
  typedef Inverted_weighted_point_2<K>   Inverted_weighted_point;
  typedef typename K::FT                 FT;
  typedef typename K::Sign               Sign;

public:
  inline Sign
  operator()(const Bitangent_line& bl,
	     const Inverted_weighted_point& v,
	     const Field_with_sqrt_tag&) const
    {
      FT a = bl.a1() + bl.a2() * CGAL::sqrt(bl.delta());
      FT b = bl.b1() + bl.b2() * CGAL::sqrt(bl.delta());
      FT c = bl.c1() + bl.c2() * CGAL::sqrt(bl.delta());
      FT r = a * v.x() + b * v.y() + c * v.p() - v.weight() * bl.d();
      return CGAL::sign(r);
    }

  inline Sign
  operator()(const Bitangent_line& bl,
	     const Inverted_weighted_point& v,
	     const Integral_domain_without_division_tag&) const
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
  inline
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2) const
  {
    if ( p1.x() == p2.x() ) {
      return p1.y() < p2.y();
    }
    return p1.x() < p2.x();
  }
};

template < class K, class MTag >
class Vertex_conflict_2
{
public:
  typedef K                                 Kernel;
  typedef MTag                              Method_tag;

  typedef typename K::Point_2               Point_2;
  typedef typename K::Site_2                Site_2;
  typedef Weighted_point_inverter_2<K>      Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>      Inverted_weighted_point;
  typedef Bitangent_line_2<K>               Bitangent_line;
  typedef Voronoi_radius_2<K>               Voronoi_radius;
  typedef typename K::FT                    FT;
  typedef typename K::Orientation           Orientation;
  typedef typename K::Sign                  Sign;
  typedef typename K::Bounded_side          Bounded_side;

  typedef Bounded_side_of_CCW_circle_2<K>   Bounded_side_of_CCW_circle;

  typedef Sign_of_distance_from_bitangent_line_2<K>
                                     Sign_of_distance_from_bitangent_line;

  typedef Sign_of_distance_from_CCW_circle_2<K>
                                         Sign_of_distance_from_CCW_circle;

private:
  inline Orientation
  orientation(const Bitangent_line& l, const Point_2& p,
	      const Field_with_sqrt_tag&) const
    {
      FT A = l.a1() * p.x() + l.b1() * p.y() + l.c1();
      FT B = l.a2() * p.x() + l.b2() * p.y() + l.c2();
      FT P = A + B * CGAL::sqrt(l.delta());
      return CGAL::sign(P);
    }

  inline Orientation
  orientation(const Bitangent_line& l, const Point_2& p,
	      const Integral_domain_without_division_tag&) const
    {
      FT A = l.a1() * p.x() + l.b1() * p.y() + l.c1();
      FT B = l.a2() * p.x() + l.b2() * p.y() + l.c2();
      return sign_a_plus_b_x_sqrt_c(A, B, l.delta());
    }

  
  inline Orientation
  orientation(const Bitangent_line& l,
	      const Inverted_weighted_point& u) const
    {
      FT A = l.a1() * u.x() / u.p() + l.b1() * u.y() / u.p() + l.c1();
      FT B = l.a2() * u.x() / u.p() + l.b2() * u.y() / u.p() + l.c2();
      FT P = A + B * CGAL::sqrt(l.delta());
      return CGAL::sign(P);
    }

public:
  typedef Sign                result_type;
  typedef Site_2              argument_type;

  inline
  Sign operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& p3, const Site_2& q) const
  {
#ifdef AG2_PROFILE_PREDICATES
    ag2_predicate_profiler::incircle_counter++;
#endif
    //
    Method_tag tag;
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

  inline
  Sign operator()(const Site_2& p1, const Site_2& p2,
		  const Site_2& q) const
  {
    Method_tag tag;
    //
    Bitangent_line bl_21(p2, p1);
    Sign s = Sign_of_distance_from_bitangent_line()(bl_21, q, tag);
    if ( s != ZERO ) { return s; }

    Bitangent_line bl1_perp = bl_21.perpendicular(p1.point());
    Bitangent_line bl2_perp = bl_21.perpendicular(p2.point());
    Orientation o1 = orientation(bl1_perp, q.point(), tag);
    Orientation o2 = orientation(bl2_perp, q.point(), tag);

    CGAL_assertion( o1 != COLLINEAR || o2 != COLLINEAR );
    if ( o1 == o2 ) { return POSITIVE; }
    return NEGATIVE;
  }  

};

//--------------------------------------------------------------------

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_INCIRCLE_C2_H

// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_AG2_INCIRCLE_C2_H
#define CGAL_AG2_INCIRCLE_C2_H

#include <CGAL/enum.h>
#include <CGAL/Number_type_traits.h>

#include <CGAL/predicates/Apollonius_graph_predicate_constructions_C2.h>

#include <CGAL/predicates/Ag2_bounded_side_of_ccw_circle_C2.h>


CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------

template< class K >
class Sign_of_distance_from_bitangent_line
{
public:
  typedef CGAL::Bitangent_line<K>           Bitangent_line;
  typedef typename K::Site_2                Site_2;
  typedef CGAL::Inverted_weighted_point<K>  Inverted_weighted_point;
  typedef typename K::FT                    FT;

public:

  inline Sign
  operator()(const Bitangent_line& bl, const Site_2& q,
	     const Sqrt_field_tag&) const
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
	     const Ring_tag&) const
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
class Sign_of_distance_from_CCW_circle
{
public:
  typedef CGAL::Bitangent_line<K>           Bitangent_line;
  typedef CGAL::Inverted_weighted_point<K>  Inverted_weighted_point;
  typedef typename K::FT                    FT;
public:

  inline Sign
  operator()(const Bitangent_line& bl,
	     const Inverted_weighted_point& v,
	     const Sqrt_field_tag&) const
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
	     const Ring_tag&) const
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
class Incircle_test
{
public:
  typedef K                                 Kernel;
  typedef MTag                              Method_tag;

  typedef typename K::Point_2               Point_2;
  typedef typename K::Site_2                Site_2;
  typedef CGAL::Weighted_point_inverter<K>  Weighted_point_inverter;
  typedef CGAL::Inverted_weighted_point<K>  Inverted_weighted_point;
  typedef CGAL::Bitangent_line<K>           Bitangent_line;
  typedef CGAL::Voronoi_radius<K>           Voronoi_radius;
  typedef typename K::FT                    FT;

  typedef CGAL::Bounded_side_of_CCW_circle<K>
                                               Bounded_side_of_CCW_circle;
  typedef CGAL::Sign_of_distance_from_bitangent_line<K>
                                     Sign_of_distance_from_bitangent_line;
  typedef CGAL::Sign_of_distance_from_CCW_circle<K>
                                         Sign_of_distance_from_CCW_circle;

private:
  inline Orientation
  orientation(const Bitangent_line& l, const Point_2& p,
	      const Sqrt_field_tag&) const
    {
      FT A = l.a1() * p.x() + l.b1() * p.y() + l.c1();
      FT B = l.a2() * p.x() + l.b2() * p.y() + l.c2();
      FT P = A + B * CGAL::sqrt(l.delta());
      return CGAL::sign(P);
    }

  inline Orientation
  orientation(const Bitangent_line& l, const Point_2& p,
	      const Ring_tag&) const
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
  struct Arity {};

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

    CGAL_assertion( o1 != COLLINEAR && o2 != COLLINEAR );
    if ( o1 == o2 ) { return POSITIVE; }
    return NEGATIVE;
  }  

};

//--------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_AG2_INCIRCLE_C2_H

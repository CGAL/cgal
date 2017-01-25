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




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_BASIC_PREDICATES_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_BASIC_PREDICATES_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>



#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Segment_Delaunay_graph_2/Sqrt_extension_2.h>



namespace CGAL {

namespace SegmentDelaunayGraph_2 {

template<class K>
struct Basic_predicates_C2
{
public:
  //-------------------------------------------------------------------
  // TYPES
  //-------------------------------------------------------------------

  typedef typename K::RT                  RT;
  typedef typename K::FT                  FT;
  typedef typename K::Point_2             Point_2;
  typedef typename K::Segment_2           Segment_2;
  typedef typename K::Site_2              Site_2;
  typedef typename K::Oriented_side       Oriented_side;
  typedef typename K::Comparison_result   Comparison_result;
  typedef typename K::Sign                Sign;
  typedef typename K::Orientation         Orientation;
  typedef typename K::Compute_scalar_product_2 Compute_scalar_product_2;
  typedef typename K::Boolean             Boolean;

  typedef CGAL::Sqrt_extension<RT,RT,Tag_true>       Sqrt_1;
  typedef CGAL::Sqrt_extension_2<RT>                         Sqrt_2;
  typedef CGAL::Sqrt_extension_2<Sqrt_1>                     Sqrt_3;
private:
    typedef typename Algebraic_structure_traits<RT>::Algebraic_category RT_Category;
    typedef typename Algebraic_structure_traits<FT>::Algebraic_category FT_Category;
public:
    typedef Boolean_tag<CGAL::is_same_or_derived<Field_with_sqrt_tag,RT_Category>::value>  RT_Has_sqrt;
    typedef Boolean_tag<CGAL::is_same_or_derived<Field_with_sqrt_tag,FT_Category>::value>  FT_Has_sqrt;


  class Line_2
  {
  private:
    RT a_, b_, c_;

  public:
    Line_2() : a_(1), b_(0), c_(0) {}
    Line_2(const RT& a, const RT& b, const RT& c)
      : a_(a), b_(b), c_(c) {}

    RT a() const { return a_; }
    RT b() const { return b_; }
    RT c() const { return c_; }


    Oriented_side oriented_side(const Point_2& p) const
    {
      Sign s = CGAL::sign(a_ * p.x() + b_ * p.y() + c_);
      if ( s == ZERO ) { return ON_ORIENTED_BOUNDARY; }
      return (s == POSITIVE) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
    }
  };

  class Homogeneous_point_2
  {
  private:
    RT hx_, hy_, hw_;

  public:
    Homogeneous_point_2() : hx_(0), hy_(0), hw_(1) {}
    Homogeneous_point_2(const RT& hx, const RT& hy, const RT& hw)
      :  hx_(hx), hy_(hy), hw_(hw)
    {
      CGAL_precondition( !(CGAL::is_zero(hw_)) );
    }

    Homogeneous_point_2(const Point_2& p)
      : hx_(p.x()), hy_(p.y()), hw_(1) {}

    Homogeneous_point_2(const Homogeneous_point_2& other)
      : hx_(other.hx_), hy_(other.hy_), hw_(other.hw_) {}

    RT hx() const { return hx_; }
    RT hy() const { return hy_; }
    RT hw() const { return hw_; }

    FT x() const { return hx_ / hw_; }
    FT y() const { return hy_ / hw_; }
  };

public:
  //-------------------------------------------------------------------
  // CONVERSIONS
  //-------------------------------------------------------------------

  static FT compute_sqrt(const FT& x, const Tag_true&)
  {
    return CGAL::sqrt( x );
  }

  static FT compute_sqrt(const FT& x, const Tag_false&)
  {
    return FT(  CGAL::sqrt( CGAL::to_double(x) )  );
  }

  static
  FT to_ft(const Sqrt_1& x)
  {
    FT sqrt_c = compute_sqrt( x.root(), FT_Has_sqrt() );
    return x.a0() + x.a1() * sqrt_c;
  }

  static
  FT to_ft(const Sqrt_3& x)
  {
    FT sqrt_e = compute_sqrt( to_ft(x.e()), FT_Has_sqrt() );
    FT sqrt_f = compute_sqrt( to_ft(x.f()), FT_Has_sqrt() );
    FT sqrt_ef = sqrt_e * sqrt_f;
    return to_ft(x.a()) + to_ft(x.b()) * sqrt_e
      + to_ft(x.c()) * sqrt_f + to_ft(x.d()) * sqrt_ef;
  }

public:    //    compute_supporting_line(q.supporting_segment(), a1, b1, c1);
    //    compute_supporting_line(r.supporting_segment(), a2, b2, c2);

  //-------------------------------------------------------------------
  // BASIC CONSTRUCTIONS
  //-------------------------------------------------------------------

#if 1
  static
  Line_2 compute_supporting_line(const Site_2& s)
  {
    RT a, b, c;
    compute_supporting_line(s, a, b, c);
    return Line_2(a, b, c);
  }

  static
  void compute_supporting_line(const Site_2& s,
			       RT& a, RT& b, RT& c)
  {
    a = s.source().y() - s.target().y();
    b = s.target().x() - s.source().x();
    c = s.source().x() * s.target().y() - s.target().x() * s.source().y();
  }
#else
  static
  Line_2 compute_supporting_line(const Segment_2& s)
  {
    RT a, b, c;
    compute_supporting_line(s, a, b, c);
    return Line_2(a, b, c);
  }

  static
  void compute_supporting_line(const Segment_2& s,
			       RT& a, RT& b, RT& c)
  {
    a = s.source().y() - s.target().y();
    b = s.target().x() - s.source().x();
    c = s.source().x() * s.target().y() - s.target().x() * s.source().y();
  }
#endif

  static
  Homogeneous_point_2
  compute_projection(const Line_2& l, const Point_2& p)
  {
    RT ab = l.a() * l.b();

    RT hx = CGAL::square(l.b()) * p.x()
      - ab * p.y() - l.a() * l.c();
    RT hy = CGAL::square(l.a()) * p.y()
      - ab * p.x() - l.b() * l.c();
    RT hw = CGAL::square(l.a()) + CGAL::square(l.b());

    return Homogeneous_point_2(hx, hy, hw);
  }


  static
  Homogeneous_point_2
  projection_on_line(const Line_2& l, const Point_2& p)
  {
    RT ab = l.a() * l.b();

    RT hx = CGAL::square(l.b()) * p.x()
      - ab * p.y() - l.a() * l.c();
    RT hy = CGAL::square(l.a()) * p.y()
      - ab * p.x() - l.b() * l.c();
    RT hw = CGAL::square(l.a()) + CGAL::square(l.b());

    return Homogeneous_point_2(hx, hy, hw);
  }


  static
  Homogeneous_point_2
  midpoint(const Point_2& p1, const Point_2& p2)
  {
    RT hx = p1.x() + p2.x();
    RT hy = p1.y() + p2.y();
    RT hw = RT(2);

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Homogeneous_point_2
  midpoint(const Homogeneous_point_2& p1,
           const Homogeneous_point_2& p2)
  {
    RT hx = p1.hx() * p2.hw() + p2.hx() * p1.hw();
    RT hy = p1.hy() * p2.hw() + p2.hy() * p1.hw();
    RT hw = RT(2) * p1.hw() * p2.hw();

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Line_2 compute_perpendicular(const Line_2& l, const Point_2& p)
  {
    RT a, b, c;
    a = -l.b();
    b = l.a();
    c = l.b() * p.x() - l.a() * p.y();
    return Line_2(a, b, c);
  }

  static
  Line_2 opposite_line(const Line_2& l)
  {
    return Line_2(-l.a(), -l.b(), -l.c());
  }

  static
  RT compute_squared_distance(const Point_2& p, const Point_2& q)
  {
    return CGAL::square(p.x() - q.x()) + CGAL::square(p.y() - q.y());
  }

  static
  std::pair<RT,RT>
  compute_squared_distance(const Point_2& p, const Line_2& l)
  {
    RT d2 = CGAL::square(l.a() * p.x() + l.b() * p.y() + l.c());
    RT n2 = CGAL::square(l.a()) + CGAL::square(l.b());
    return std::pair<RT,RT>(d2, n2);
  }

public:
  //-------------------------------------------------------------------
  // BASIC PREDICATES
  //-------------------------------------------------------------------
  static
  Comparison_result
  compare_squared_distances_to_line(const Line_2& l, const Point_2& p,
                                    const Point_2& q)
  {
    RT d2_lp = CGAL::square(l.a() * p.x() + l.b() * p.y() + l.c());
    RT d2_lq = CGAL::square(l.a() * q.x() + l.b() * q.y() + l.c());

    return CGAL::compare(d2_lp, d2_lq);
  }

  static
  Comparison_result
  compare_squared_distances_to_lines(const Point_2& p,
				     const Line_2& l1,
                                     const Line_2& l2)
  {
    RT d2_l1 = CGAL::square(l1.a() * p.x() + l1.b() * p.y() + l1.c());

    RT d2_l2 = CGAL::square(l2.a() * p.x() + l2.b() * p.y() + l2.c());

    RT n1 = CGAL::square(l1.a()) + CGAL::square(l1.b());
    RT n2 = CGAL::square(l2.a()) + CGAL::square(l2.b());

    return CGAL::compare(d2_l1 * n2, d2_l2 * n1);
  }

  static
  Oriented_side
  oriented_side_of_line(const Line_2& l, const Point_2& p)
  {
    return CGAL::sign(l.a() * p.x() + l.b() * p.y() + l.c());
  }

  static
  Oriented_side
  oriented_side_of_line(const Line_2& l, const Homogeneous_point_2& p)
  {
    Sign s1 =
      CGAL::sign(l.a() * p.hx() + l.b() * p.hy() + l.c() * p.hw());
    Sign s_hw = CGAL::sign(p.hw());

    return s1 * s_hw;
  }


  static
  bool is_on_positive_halfspace(const Line_2& l, const Segment_2& s)
  {
    Oriented_side os1, os2;

    os1 = oriented_side_of_line(l, s.source());
    os2 = oriented_side_of_line(l, s.target());

    return ( (os1 == ON_POSITIVE_SIDE && os2 != ON_NEGATIVE_SIDE) ||
	     (os1 != ON_NEGATIVE_SIDE && os2 == ON_POSITIVE_SIDE) );
  }

};


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_BASIC_PREDICATES_C2_H

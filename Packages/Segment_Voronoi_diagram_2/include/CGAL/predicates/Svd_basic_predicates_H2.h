// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/predicates/Svd_basic_predicates_H2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_BASIC_PREDICATES_H2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_BASIC_PREDICATES_H2_H


#include <CGAL/enum.h>
#include <CGAL/predicates/Square_root_1.h>
#include <CGAL/predicates/Square_root_2.h>



CGAL_BEGIN_NAMESPACE


template<class K>
struct Svd_basic_predicates_H2
{
public:
  //-------------------------------------------------------------------
  // TYPES
  //-------------------------------------------------------------------

  typedef typename K::RT         RT;
  typedef typename K::FT         FT;
  typedef typename K::Point_2    Point_2;
  typedef typename K::Segment_2  Segment_2;
  typedef typename K::Site_2     Site_2;

  typedef CGAL::Square_root_1<RT>       Sqrt_1;
  typedef CGAL::Square_root_2<RT>       Sqrt_2;
  typedef CGAL::Square_root_2<Sqrt_1>   Sqrt_3;

  typedef Point_2                Homogeneous_point_2;

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

  };

#if 0
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
#endif

public:
  //-------------------------------------------------------------------
  // CONVERSIONS
  //-------------------------------------------------------------------

  static
  FT to_ft(const Sqrt_1& x)
  {
    FT sqrt_c = CGAL::sqrt(x.c());
    return x.a() + x.b() * sqrt_c;
  }

  static
  FT to_ft(const Sqrt_3& x)
  {
    FT sqrt_e = CGAL::sqrt(x.e());
    FT sqrt_f = CGAL::sqrt(x.f());
    FT sqrt_ef = sqrt_e * sqrt_f;
    return x.a() + x.b() * sqrt_e +  x.c() * sqrt_f + x.d() * sqrt_ef;
  }

public:
  //-------------------------------------------------------------------
  // BASIC CONSTRUCTIONS
  //-------------------------------------------------------------------

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
    RT x1 = s.source().hx();
    RT y1 = s.source().hy();
    RT w1 = s.source().hw();

    if ( w1 < 0 ) {
      x1 = -x1;
      y1 = -y1;
      w1 = -w1;
    }

    RT x2 = s.target().hx();
    RT y2 = s.target().hy();
    RT w2 = s.target().hw();

    if ( w2 < 0 ) {
      x2 = -x2;
      y2 = -y2;
      w2 = -w2;
    }

    a = y1 * w2 - y2 * w1;
    b = x2 * w1 - x1 * w2;
    c = x1 * y2 - x2 * y1;
  }

  static
  Homogeneous_point_2
  compute_projection(const Line_2& l, const Point_2& p)
  {
    Homogeneous_point_2 hp;

    RT ab = l.a() * l.b();

    RT hx = CGAL::square(l.b()) * p.hx()
      - ab * p.hy() - l.a() * l.c() * p.hw();
    RT hy = CGAL::square(l.a()) * p.hy()
      - ab * p.hx() - l.b() * l.c() * p.hw();
    RT hw = CGAL::square(l.a()) + CGAL::square(l.b());

    return Homogeneous_point_2(hx, hy, hw * p.hw());
  }

  static
  Homogeneous_point_2
  projection_on_line(const Line_2& l, const Point_2& p)
  {
    return compute_projection(l, p);
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
    RT hx = p.hx();
    RT hy = p.hy();
    RT hw = p.hw();

    if ( w < 0 ) {
      hx = -hx;
      hy = -hy;
      hw = -hw;
    }

    RT a, b, c;
    a = -l.b() * hw;
    b = l.a() * hw;
    c = l.b() * hx - l.a() * hy;
    return Line_2(a, b, c);
  }

  static
  Line_2 opposite_line(const Line_2& l)
  {
    return Line_2(-l.a(), -l.b(), -l.c());
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
    RT d2_lp =
      CGAL::square(l.a() * p.hx() + l.b() * p.hy() + l.c() * p.hw());
    RT d2_lq =
      CGAL::square(l.a() * q.hx() + l.b() * q.hy() + l.c() * p.hw());

    RT p_hw_sq = CGAL::square( p.hw() );
    RT q_hw_sq = CGAL::square( q.hw() );

    return CGAL::compare(d2_lp * q_hw_sq, d2_lq * p_hw_sq);
  }

  static
  Comparison_result
  compare_squared_distances_to_lines(const Point_2& p,
				     const Line_2& l1,
                                     const Line_2& l2)
  {
    RT hx = p.hx();
    RT hy = p.hy();
    RT hw = p.hw();

    RT d2_l1 = CGAL::square(l1.a() * hx + l1.b() * hy + l1.c() * hw);

    RT d2_l2 = CGAL::square(l2.a() * hx + l2.b() * hy + l2.c() * hw);

    RT n1 = CGAL::square(l1.a()) + CGAL::square(l1.b());
    RT n2 = CGAL::square(l2.a()) + CGAL::square(l2.b());

    return CGAL::compare(d2_l1 * n2, d2_l2 * n1);
  }

  static
  Oriented_side
  oriented_side_of_line(const Line_2& l, const Homogeneous_point_2& p)
  {
    Sign s1 =
      CGAL::sign(l.a() * p.hx() + l.b() * p.hy() + l.c() * p.hw());
    Sign s_hw = CGAL::sign(p.hw());

    Sign s = Sign(s1 * s_hw);

    if ( s == ZERO ) { return ON_ORIENTED_BOUNDARY; }
    return ( s == POSITIVE ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
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


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_BASIC_PREDICATES_C2_H

#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_BASIC_PREDICATES_C2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_BASIC_PREDICATES_C2_H


#include <CGAL/enum.h>
#include <CGAL/predicates/Square_root_1.h>
#include <CGAL/predicates/Square_root_2.h>



CGAL_BEGIN_NAMESPACE


template<class K>
struct Svd_basic_predicates_C2
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

  class Homogeneous_point_2
  {
  private:
    RT hx_, hy_, hw_;

  public:
    Homogeneous_point_2() : hx_(0), hy_(0), hw_(1) {}
    Homogeneous_point_2(const RT& hx, const RT& hy, const RT& hw)
      :  hx_(hx), hy_(hy), hw_(hw)
    {
      CGAL_precondition( !(CGAL_NTS is_zero(hw_)) );
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

  static
  FT to_ft(const Sqrt_1& x)
  {
    FT sqrt_c = CGAL_NTS sqrt(x.c());
    return x.a() + x.b() * sqrt_c;
  }

  static
  FT to_ft(const Sqrt_3& x)
  {
    FT sqrt_e = CGAL_NTS sqrt(x.e());
    FT sqrt_f = CGAL_NTS sqrt(x.f());
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
    RT a = s.source().y() - s.target().y();
    RT b = s.target().x() - s.source().x();
    RT c = s.source().x() * s.target().y()
      - s.target().x() * s.source().y();
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

  static
  Homogeneous_point_2
  compute_projection(const Line_2& l, const Point_2& p)
  {
    Homogeneous_point_2 hp;

    RT ab = l.a() * l.b();

    RT hx = CGAL_NTS square(l.b()) * p.x()
      - ab * p.y() - l.a() * l.c();
    RT hy = CGAL_NTS square(l.a()) * p.y()
      - ab * p.x() - l.b() * l.c();
    RT hw = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());

    return Homogeneous_point_2(hx, hy, hw);
  }


  static
  Homogeneous_point_2
  projection_on_line(const Line_2& l, const Point_2& p)
  {
    RT ab = l.a() * l.b();

    RT hx = CGAL_NTS square(l.b()) * p.x()
      - ab * p.y() - l.a() * l.c();
    RT hy = CGAL_NTS square(l.a()) * p.y()
      - ab * p.x() - l.b() * l.c();
    RT hw = CGAL_NTS square(l.a()) + CGAL_NTS square(l.b());

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

public:
  //-------------------------------------------------------------------
  // BASIC PREDICATES
  //-------------------------------------------------------------------
  static
  Comparison_result
  compare_squared_distances_to_line(const Line_2& l, const Point_2& p,
                                    const Point_2& q)
  {
    RT d2_lp = CGAL_NTS square(l.a() * p.x() + l.b() * p.y() + l.c());
    RT d2_lq = CGAL_NTS square(l.a() * q.x() + l.b() * q.y() + l.c());

    return CGAL_NTS compare(d2_lp, d2_lq);
  }

  static
  Comparison_result
  compare_squared_distances_to_lines(const Point_2& p,
				     const Line_2& l1,
                                     const Line_2& l2)
  {
    RT d2_l1 =
      CGAL_NTS square(l1.a() * p.x() + l1.b() * p.y() + l1.c());

    RT d2_l2 =
      CGAL_NTS square(l2.a() * p.x() + l2.b() * p.y() + l2.c());

    RT n1 = CGAL_NTS square(l1.a()) + CGAL_NTS square(l1.b());
    RT n2 = CGAL_NTS square(l2.a()) + CGAL_NTS square(l2.b());

    return CGAL_NTS compare(d2_l1 * n2, d2_l2 * n1);
  }

  static
  Oriented_side
  oriented_side_of_line(const Line_2& l, const Point_2& p)
  {
    Sign s = CGAL_NTS sign(l.a() * p.x() + l.b() * p.y() + l.c());
    if ( s == ZERO ) { return ON_ORIENTED_BOUNDARY; }
    return ( s == POSITIVE ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  }

  static
  Oriented_side
  oriented_side_of_line(const Line_2& l, const Homogeneous_point_2& p)
  {
    Sign s1 =
      CGAL_NTS sign(l.a() * p.hx() + l.b() * p.hy() + l.c() * p.hw());
    Sign s_hw = CGAL_NTS sign(p.hw());

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


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_BASIC_PREDICATES_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_BASIC_PREDICATES_C2_H


#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Sqrt_extension_2.h>
#include <CGAL/Orientation_Linf_2.h>

#include <CGAL/Polychain_2.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

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
  typedef typename K::Bounded_side        Bounded_side;
  typedef typename K::Comparison_result   Comparison_result;
  typedef typename K::Sign                Sign;
  typedef typename K::Orientation         Orientation;
  typedef typename CGAL::OrientationLinf  OrientationLinf;
  typedef typename K::Compute_scalar_product_2 
                        Compute_scalar_product_2;
  typedef typename K::Boolean             Boolean;
  typedef typename K::Direction_2         Direction_2;
  typedef typename K::Vector_2            Vector_2;
  typedef typename K::Compare_x_2         Compare_x_2;
  typedef typename K::Compare_y_2         Compare_y_2;

  typedef typename CGAL::Polychainline_2<K> Polychainline_2;   

  typedef CGAL::Sqrt_extension<RT,RT,Tag_true>     Sqrt_1;
  typedef CGAL::Sqrt_extension_2<RT>       Sqrt_2;
  typedef CGAL::Sqrt_extension_2<Sqrt_1>   Sqrt_3;
private:
    typedef typename Algebraic_structure_traits<RT>::Algebraic_category RT_Category;
    typedef typename Algebraic_structure_traits<FT>::Algebraic_category FT_Category;
public:
    typedef Boolean_tag<CGAL::is_same_or_derived<Field_with_sqrt_tag,RT_Category>::value>  RT_Has_sqrt;
    typedef Boolean_tag<CGAL::is_same_or_derived<Field_with_sqrt_tag,FT_Category>::value>  FT_Has_sqrt;
 
  static const RT_Has_sqrt& rt_has_sqrt() {
    static RT_Has_sqrt has_sqrt;
    return has_sqrt;
  }

  static const FT_Has_sqrt& ft_has_sqrt() {
    static FT_Has_sqrt has_sqrt;
    return has_sqrt;
  }
				   

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
    FT sqrt_c = compute_sqrt( x.root(), ft_has_sqrt() );
    return x.a0() + x.a1() * sqrt_c;
  }

  static
  FT to_ft(const Sqrt_3& x)
  {
    FT sqrt_e = compute_sqrt( to_ft(x.e()), ft_has_sqrt() );
    FT sqrt_f = compute_sqrt( to_ft(x.f()), ft_has_sqrt() );
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

  // compute horizontal line l that goes through p,
  // and leaves q on the oriented side s
  // s: has to be either +1 or -1 (not 0)
  // q: should not be on line l
  static 
  Line_2
  compute_horizontal_side_line(
      const Point_2& p, const Point_2& q, Oriented_side s)  
  {
    CGAL_precondition(s != ON_ORIENTED_BOUNDARY);
    
    RT b, c;

    b = RT(1);
    c = - p.y();

    // direction is (1, 0)

    Compare_y_2 cmpy;
    if (((cmpy(q, p) == LARGER) and
         (s == ON_NEGATIVE_SIDE)   ) or
        ((cmpy(q, p) == SMALLER) and
         (s == ON_POSITIVE_SIDE)   )   ) {
      b = -b;
      c = -c;
    }
    return Line_2(RT(0), b, c);
  }

  // compute vertical line l that goes through p,
  // and leaves q on the oriented side s
  // s: has to be either +1 or -1 (not 0)
  // q: should not be on line l
  static 
  Line_2
  compute_vertical_side_line(
      const Point_2& p, const Point_2& q, Oriented_side s)  
  {
    CGAL_precondition(s != ON_ORIENTED_BOUNDARY);
    
    RT a, c;

    a = RT(1);
    c = - p.x();

    // direction is (0, -1)

    Compare_x_2 cmpx;
    if (((cmpx(q, p) == LARGER) and
         (s == ON_NEGATIVE_SIDE)   ) or
        ((cmpx(q, p) == SMALLER) and
         (s == ON_POSITIVE_SIDE)   )   ) {
      a = -a;
      c = -c;
    }
    return Line_2(a, RT(0), c);
  }


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
  compute_linf_projection_hom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( not l.is_degenerate() );
    CGAL_precondition( 
        (CGAL::sign(l.a()) != ZERO) or (CGAL::sign(l.b()) != ZERO) );

    Sign signa = CGAL::sign(l.a()); 
    Sign signb = CGAL::sign(l.b());

    RT hx, hy, hw;

    if (signa == ZERO) { // l is horizontal
      // l.a() == 0  =>  y = -c/b
      hx = p.x() * l.b();
      hy = - l.c();
      hw = l.b();
    } else if (signb == ZERO) { // l is vertical
      // l.b() == 0  =>  x = -c/a
      hx = - l.c();
      hy = p.y() * l.a();
      hw = l.a();
    } else {
      // here both l.a() and l.b() are non-zero
      if ( signa == signb ) {
        hx = l.b() * ( p.x() - p.y() ) - l.c();
        hy = l.a() * ( p.y() - p.x() ) - l.c();
        hw = l.a() + l.b();
      } else { // signa != signb
        hx = -l.b() * ( p.x() + p.y() ) - l.c();
        hy = l.a() * ( p.x() + p.y() ) + l.c();
        hw = l.a() - l.b();
      }
    }

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Point_2
  compute_linf_projection_nonhom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( not l.is_degenerate() );
    CGAL_precondition( 
        (CGAL::sign(l.a()) != ZERO) or (CGAL::sign(l.b()) != ZERO) );

    Sign signa = CGAL::sign(l.a()); 
    Sign signb = CGAL::sign(l.b());

    RT hx, hy, hw;

    if (signa == ZERO) { // l is horizontal
      // l.a() == 0  =>  y = -c/b
      hx = p.x() * l.b();
      hy = - l.c();
      hw = l.b();
    } else if (signb == ZERO) { // l is vertical
      // l.b() == 0  =>  x = -c/a
      hx = - l.c();
      hy = p.y() * l.a();
      hw = l.a();
    } else {
      // here both l.a() and l.b() are non-zero
      if ( signa == signb ) {
        hx = l.b() * ( p.x() - p.y() ) - l.c();
        hy = l.a() * ( p.y() - p.x() ) - l.c();
        hw = l.a() + l.b();
      } else { // signa != signb
        hx = -l.b() * ( p.x() + p.y() ) - l.c();
        hy = l.a() * ( p.x() + p.y() ) + l.c();
        hw = l.a() - l.b();
      }
    }

    return Point_2(hx, hy, hw);
  }

  static
  Homogeneous_point_2
  compute_horizontal_projection_hom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( not l.is_horizontal() );
    CGAL_precondition( 
        (CGAL::sign(l.a()) != ZERO) );

    RT hx, hy, hw;

    hx = - l.b() * p.y() - l.c();
    hy = p.y() * l.a();
    hw = l.a();

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Point_2
  compute_horizontal_projection(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( not l.is_horizontal() );
    CGAL_precondition( 
        (CGAL::sign(l.a()) != ZERO) );

    RT hx, hy, hw;

    hx = - l.b() * p.y() - l.c();
    hy = p.y() * l.a();
    hw = l.a();

    return Point_2(hx, hy, hw);
  }
  
  static
  Homogeneous_point_2
  compute_vertical_projection_hom(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( not l.is_horizontal() );
    CGAL_precondition( 
        (CGAL::sign(l.b()) != ZERO) );

    RT hx, hy, hw;

    hx = p.x() * l.b(); 
    hy = - l.a() * p.x() - l.c();
    hw = l.b();

    return Homogeneous_point_2(hx, hy, hw);
  }

  static
  Point_2
  compute_vertical_projection(const Line_2& l, const Point_2& p)
  {
    //CGAL_precondition( not l.is_horizontal() );
    CGAL_precondition( 
        (CGAL::sign(l.b()) != ZERO) );

    RT hx, hy, hw;

    hx = p.x() * l.b(); 
    hy = - l.a() * p.x() - l.c();
    hw = l.b();

    return Point_2(hx, hy, hw);
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

  // compute_cw_perpendicular is opposite of compute_perpendicular
  static
  Line_2 compute_cw_perpendicular(const Line_2& l, const Point_2& p)
  {
    RT a, b, c;
    a = l.b();
    b = -l.a();
    c = -l.b() * p.x() + l.a() * p.y();
    return Line_2(a, b, c);
  }

  static
  Line_2 compute_linf_perpendicular(const Line_2& l, const Point_2& p)
  {
    RT a, b, c;
    a = RT( - CGAL::sign(l.b()) );
    b = RT( CGAL::sign(l.a()) );
    c = - a * p.x() - b * p.y();
    return Line_2(a, b, c);
  }

  static
  Line_2 opposite_line(const Line_2& l)
  {
    return Line_2(-l.a(), -l.b(), -l.c());
  }

  // philaris: similar to compute_supporting_line
  static 
  Line_2 compute_line_from_to(const Point_2& p, const Point_2&q)
  {
    RT a, b, c;
    a = p.y() - q.y();
    b = q.x() - p.x();

    CGAL_assertion((CGAL::sign(a) != ZERO) or 
                   (CGAL::sign(b) != ZERO))   ;

    c = p.x() * q.y() - q.x() * p.y();

    return Line_2(a, b, c);
  }

  static
  RT compute_linf_distance(const Point_2& p, const Point_2& q)
  {
    return CGAL::max(
              CGAL::abs(p.x() - q.x()), 
              CGAL::abs(p.y() - q.y()));
  }

  static
  void compute_intersection_of_lines(
      const Line_2& l1, const Line_2& l2,
      RT& hx, RT& hy, RT& hw)
  {
    hx = l1.b() * l2.c() - l1.c() * l2.b();
    hy = l1.c() * l2.a() - l1.a() * l2.c();
    hw = l1.a() * l2.b() - l1.b() * l2.a();
  }

public:
  //-------------------------------------------------------------------
  // BASIC PREDICATES
  //-------------------------------------------------------------------
  static
  Comparison_result
  compare_linf_distances_to_line(const Line_2& l, const Point_2& p,
                                    const Point_2& q)
  {
    Homogeneous_point_2 hp = compute_linf_projection_hom(l, p);
    Homogeneous_point_2 hq = compute_linf_projection_hom(l, q);

    RT dlp = CGAL::max(CGAL::abs(p.x() - hp.x()),
                       CGAL::abs(p.y() - hp.y()));

    RT dlq = CGAL::max(CGAL::abs(q.x() - hq.x()),
                       CGAL::abs(q.y() - hq.y()));

    Comparison_result crude = CGAL::compare(dlp, dlq);

    if (crude != EQUAL) {
      return crude;
    } else {
      std::cout << "compare_linf_distances_to_line refining" 
        << std::endl;
      return crude;
    }
  }

  static
  Comparison_result
  compare_linf_distances_to_lines(const Point_2& p,
				     const Line_2& l1,
                                     const Line_2& l2)
  {
    Homogeneous_point_2 hl1 = compute_linf_projection_hom(l1, p);
    Homogeneous_point_2 hl2 = compute_linf_projection_hom(l2, p);

    RT dl1p = CGAL::max(CGAL::abs(hl1.x() - p.x()),
                        CGAL::abs(hl1.y() - p.y()));

    RT dl2p = CGAL::max(CGAL::abs(hl2.x() - p.x()),
                        CGAL::abs(hl2.y() - p.y()));

    Comparison_result crude = CGAL::compare(dl1p, dl2p);

    if (crude != EQUAL) {
      return crude;
    } else {
      std::cout << "compare_linf_distances_to_lines refining" 
        << std::endl;
      return crude;
    }
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

  static 
  Comparison_result
  compare_distance_to_point_linf(
      const Point_2& p, const Point_2& q, const Point_2& r)
  {
    Comparison_result retval;
    retval = 
      CGAL::compare(
        CGAL::max( CGAL::abs(p.x()-q.x()), CGAL::abs(p.y()-q.y()) ),
        CGAL::max( CGAL::abs(p.x()-r.x()), CGAL::abs(p.y()-r.y()) )
        );
    if (retval == EQUAL) {
      std::cout << "debug cmpdistlinf maybe break ties" << std::endl;
      // here, p might be on the 2-dimensional bisector(q,r),
      // therefore we have to break ties, based on one coordinate
      if (CGAL::compare(q.x(), r.x()) == EQUAL) {
        std::cout << "debug cmpdistlinf try breaking with y" << std::endl;
        return CGAL::compare( CGAL::abs(p.y()-q.y()),
                              CGAL::abs(p.y()-r.y()) ) ;
      } else if (CGAL::compare(q.y(), r.y()) == EQUAL) {
        std::cout << "debug cmpdistlinf try breaking with x" << std::endl;
        return CGAL::compare( CGAL::abs(p.x()-q.x()),
                              CGAL::abs(p.x()-r.x()) ) ;
      }
    }
    return retval;
  }

  static 
  Bounded_side
  bounded_side_of_bbox(
      const Point_2& p, const Point_2& q, const Point_2& r)
  {
    // precondition: p, q, r should be monotone points.
    // Predicate bounded_side_of_bbox (P_bbox) returns:
    //  0 if p = q,
    //  0 if r = p or r = q (ON_BOUNDARY),
    // -1 if r is strictly outside the bounding box of p,q   
    //    (ON_UNBOUNDED_SIDE),
    // +1 if r is strictly inside the bounding box of p,q
    //    (ON_BOUNDED_SIDE).
    // If p and q are on the same vertical or horizontal 
    // line but are not the same point, then the bounding
    // box of p and q degenerates to the line segment pq.
    
    std::cout << "debug bounded_side_of_bbox (p q r) = ("
      << p << ") (" << q << ") (" << r << ")" << std::endl; 

    if ((CGAL::compare(p.x(), q.x()) == EQUAL) and
        (CGAL::compare(p.y(), q.y()) == EQUAL)    ) {
      return ON_BOUNDARY;
    }

    Comparison_result cmpxpr, cmpxrq, cmpypr, cmpyrq;

    cmpxpr = CGAL::compare(p.x(), r.x());
    cmpxrq = CGAL::compare(r.x(), q.x());
    cmpypr = CGAL::compare(p.y(), r.y());
    cmpyrq = CGAL::compare(r.y(), q.y());

    Comparison_result comp =
      CGAL::compare(cmpxpr*cmpxrq + cmpypr*cmpyrq, 0);

    std::cout << "debug bounded_side_of_bbox returns ";

    switch(comp) {
      case SMALLER:
        std::cout << "ON_UNBOUNDED_SIDE" << std::endl;
        return ON_UNBOUNDED_SIDE;
      case EQUAL:
        std::cout << "ON_BOUNDARY" << std::endl;
        return ON_BOUNDARY;
      case LARGER:
        std::cout << "ON_BOUNDED_SIDE" << std::endl;
        return ON_BOUNDED_SIDE;
      default:
        std::cout << "error: should never reach here";
        CGAL_assertion( false );
        return ON_BOUNDARY;
    }
  }


  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the positive halfplane of oriented line l
  static 
  Boolean
  intersects_segment_positive_halfplane(const Site_2 & s, 
      const Line_2 & l)
  {
    Segment_2 seg = s.segment();
    
    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    std::cout << "debug: intersects_segment_positive_halfplane "
      << "s=" << s 
      << " l=" << l.a() << " " << l.b() << " " << l.c()
      << std::endl;

    Oriented_side oslsrc = oriented_side_of_line(l, ssrc);
    Oriented_side osltrg = oriented_side_of_line(l, strg);

    std::cout << "debug: intersects_segment_positive_halfplane "
      << "oslsrc=" << oslsrc << " osltrg=" << osltrg 
      << std::endl;

    if ((oslsrc == ON_POSITIVE_SIDE) or
        (osltrg == ON_POSITIVE_SIDE)   )
    {
      return true;
    } else {
      return false;
    }
  }


  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the negative halfplane of oriented line l
  static 
  Boolean
  intersects_segment_negative_halfplane(const Site_2 & s, 
      const Line_2 & l)
  {
    Segment_2 seg = s.segment();
    
    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    std::cout << "debug: intersects_segment_negative_halfplane "
      << "s=" << s 
      << " l=" << l.a() << " " << l.b() << " " << l.c()
      << std::endl;

    Oriented_side oslsrc = oriented_side_of_line(l, ssrc);
    Oriented_side osltrg = oriented_side_of_line(l, strg);

    std::cout << "debug: intersects_segment_negative_halfplane "
      << "oslsrc=" << oslsrc << " osltrg=" << osltrg 
      << std::endl;

    if ((oslsrc == ON_NEGATIVE_SIDE) or
        (osltrg == ON_NEGATIVE_SIDE)   )
    {
      return true;
    } else {
      return false;
    }
  }


  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // the only finite corner of the infinite box is corner
  // and if you traverse the infinite box ccw, then
  // you meet points in that order: q, corner, p
  static 
  Boolean
  intersects_segment_interior_inf_box(const Site_2 & s, 
      const Point_2 & q, const Point_2 & corner, 
      const Point_2 & p)
  {
    CGAL_assertion(s.is_segment());
    Segment_2 seg = s.segment();

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    Line_2 lqc = compute_line_from_to(q, corner);
    Line_2 lcp = compute_line_from_to(corner, p);

    Oriented_side os_lqc_ssrc = oriented_side_of_line(lqc, ssrc);
    Oriented_side os_lcp_ssrc = oriented_side_of_line(lcp, ssrc);

    Oriented_side os_lqc_strg = oriented_side_of_line(lqc, strg);
    Oriented_side os_lcp_strg = oriented_side_of_line(lcp, strg);

    std::cout << "debug qcp= (" << q << ") (" << corner 
      << ") (" << p << ")" << std::endl; 

    if (((os_lqc_ssrc == ON_POSITIVE_SIDE) and 
         (os_lcp_ssrc == ON_POSITIVE_SIDE)) or
        ((os_lcp_strg == ON_POSITIVE_SIDE) and 
         (os_lqc_strg == ON_POSITIVE_SIDE))   ) {
      std::cout << "debug is_segment_inside_inf_box " 
        << "endpoint inside" << std::endl;
      return true;
    } else {
      // here you have to check if the interior is inside

      std::cout << "debug is_segment_inside_inf_box " 
        << "try for interior to be inside" << std::endl;

      // in fact, here you can intersect the segment
      // with the ray starting from corner and going to the
      // direction of the center of the infinite box

      Compare_x_2 cmpx;
      Compare_y_2 cmpy;

      Comparison_result cmpxpq = cmpx(p,q);
      Comparison_result cmpypq = cmpy(p,q);

      RT one(1);

      Point_2 displaced ( corner.x() + (-cmpypq)*one ,
                          corner.y() + cmpxpq * one   );

      Line_2 l = compute_line_from_to(corner, displaced);

      Line_2 lseg = compute_supporting_line(s);

      RT hx, hy, hw;

      compute_intersection_of_lines(l, lseg, hx, hy, hw);

      if (CGAL::sign(hw) == ZERO) {
        return false;
      } else {
        Point_2 ip ( hx/hw, hy/hw);
        Oriented_side os_lqc_ip = oriented_side_of_line(lqc, ip);
        Oriented_side os_lcp_ip = oriented_side_of_line(lcp, ip);

        Compare_x_2 cmpx;
        Compare_y_2 cmpy;

        Comparison_result cmpxsrcip = cmpx(ssrc, ip);
        Comparison_result cmpysrcip = cmpy(ssrc, ip);
        Comparison_result cmpxiptrg = cmpx(ip, strg);
        Comparison_result cmpyiptrg = cmpy(ip, strg);

        // philaris: to check
        Boolean is_ip_inside_segment = 
          (CGAL::sign(cmpxsrcip * cmpxiptrg + 
                      cmpysrcip * cmpyiptrg   )) == POSITIVE;
        

        if ((os_lqc_ip == ON_POSITIVE_SIDE) and 
            (os_lcp_ip == ON_POSITIVE_SIDE) and 
            is_ip_inside_segment ) {
          return true;
        } else {
          return false;
        }
      }
    }
  } // end of intersects_segment_interior_inf_box



  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // this infinite box is a 90 degree wedge defined
  // by the intersection of the halfplanes
  // with supporting lines lhor and lver, where
  // the halfplanes are both on the positive or negative
  // sides of the supporting lines
  static 
  Boolean
  intersects_segment_side_of_wedge(const Site_2 & s, 
      const Line_2 & lhor, const Line_2 & lver,
      Oriented_side orside) 
  {
    CGAL_assertion(s.is_segment());
    Segment_2 seg = s.segment();

    std::cout << "debug sofw s=" << s 
      << " orside=" << orside << std::endl;

    Point_2 ssrc = seg.source();
    Point_2 strg = seg.target();

    Oriented_side os_lhor_ssrc = oriented_side_of_line(lhor, ssrc);
    Oriented_side os_lver_ssrc = oriented_side_of_line(lver, ssrc);

    Oriented_side os_lhor_strg = oriented_side_of_line(lhor, strg);
    Oriented_side os_lver_strg = oriented_side_of_line(lver, strg);

    if (((os_lhor_ssrc == orside) and 
         (os_lver_ssrc == orside)) or
        ((os_lhor_strg == orside) and 
         (os_lver_strg == orside))   ) {
      std::cout << "debug intersects_segment_side_of_wedge " 
        << "endpoint inside" << std::endl;
      return true;
    } else {
      // here we have to check if the interior is inside

      std::cout << "debug intersects_segment_side_of_wedge " 
        << "try for interior to be inside" << std::endl;

      // in fact, here you can intersect the segment
      // with the ray starting from corner and going to the
      // direction of the center of the infinite box

      // corner has homogenuous coordinates cx, cy, cw
      RT cx, cy, cw;
      compute_intersection_of_lines(lhor, lver, cx, cy, cw);

      CGAL_assertion( CGAL::sign(cw) != ZERO );

      Point_2 corner ( cx, cy, cw );

      std::cout << "debug corner=" << corner << std::endl;
 
      RT one(1);

      Point_2 displaced ( 
          corner.x() + ( (+orside)*CGAL::sign(lver.a()) ) * one ,
          corner.y() +   (+orside)*CGAL::sign(lhor.b())   * one   );

      std::cout << "debug displaced=" << displaced << std::endl;

      Line_2 l = compute_line_from_to(corner, displaced);

      Line_2 lseg = compute_supporting_line(s);

      RT hx, hy, hw;

      std::cout << "debug: intersects_segment_side_of_wedge "
        << " l=" << l.a() << " " << l.b() << " " << l.c()
        << " lseg=" << lseg.a() << " " << lseg.b() << " " << lseg.c()
        << std::endl;

      compute_intersection_of_lines(l, lseg, hx, hy, hw);

      if (CGAL::sign(hw) == ZERO) {
        std::cout << "debug l and lseg are parallel" << std::endl;
        return false;
      } else {
        Point_2 ip ( hx, hy, hw );
        std::cout << "debug ip=" << ip << std::endl;
        Oriented_side os_lhor_ip = oriented_side_of_line(lhor, ip);
        Oriented_side os_lver_ip = oriented_side_of_line(lver, ip);

        Compare_x_2 cmpx;
        Compare_y_2 cmpy;

        Comparison_result cmpxsrcip = cmpx(ssrc, ip);
        Comparison_result cmpysrcip = cmpy(ssrc, ip);
        Comparison_result cmpxiptrg = cmpx(ip, strg);
        Comparison_result cmpyiptrg = cmpy(ip, strg);

        // philaris: to check
        Boolean is_ip_inside_segment = 
          (CGAL::sign(cmpxsrcip * cmpxiptrg + 
                      cmpysrcip * cmpyiptrg   )) == POSITIVE;

        if ((os_lhor_ip == orside) and
            (os_lver_ip == orside) and 
            is_ip_inside_segment      ) {
          return true;
        } else {
          return false;
        }
      }
    }
  } // end of intersects_segment_side_of_wedge

  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // this infinite box is a 90 degree wedge defined
  // by the intersection of the positive halfplanes
  // with supporting lines lhor and lver
  static 
  Boolean
  intersects_segment_positive_of_wedge(const Site_2 & s, 
      const Line_2 & lhor, const Line_2 & lver) 
  {
    return intersects_segment_side_of_wedge(
        s, lhor, lver, ON_POSITIVE_SIDE);
  }

  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the interior of the following infinite box:
  // this infinite box is a 90 degree wedge defined
  // by the intersection of the positive halfplanes
  // with supporting lines lhor and lver
  static 
  Boolean
  intersects_segment_negative_of_wedge(const Site_2 & s, 
      const Line_2 & lhor, const Line_2 & lver) 
  {
    return intersects_segment_side_of_wedge(
        s, lhor, lver, ON_NEGATIVE_SIDE);
  }


  // returns true if and only if 
  // the interior of s has non-empty intersection
  // with the interior of the bounding box of q, p
  // precondition: the bounding box should be non-trivial,
  // i.e., it should not be a segment
  static 
  Boolean
  intersects_segment_interior_bbox(const Site_2 & s, 
      const Point_2 & q, 
      const Point_2 & p)
  {
    Compare_x_2 cmpx;
    Compare_y_2 cmpy;
    CGAL_precondition(cmpx(p,q) != EQUAL);
    CGAL_precondition(cmpy(p,q) != EQUAL);
    CGAL_precondition(s.is_segment());

    Point_2 corner1 ( p.x(), q.y());
    Point_2 corner2 ( q.x(), p.y());

    if (CGAL::orientation( q, corner1, p ) == LEFT_TURN) {
      return intersects_segment_interior_inf_box(s, q, corner1, p)
         and intersects_segment_interior_inf_box(s, p, corner2, q);
    } else {
      return intersects_segment_interior_inf_box(s, q, corner2, p)
         and intersects_segment_interior_inf_box(s, p, corner1, q);
    }
  } // end of intersects_segment_interior_bbox
};



} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_BASIC_PREDICATES_C2_H

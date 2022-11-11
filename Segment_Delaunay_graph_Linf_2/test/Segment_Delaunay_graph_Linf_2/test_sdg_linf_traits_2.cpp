
#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_algebraic.h>

#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>




typedef CGAL::Simple_cartesian<double>      Double_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<true> >   IT_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >  IF_Kernel;
typedef CGAL::Simple_cartesian<double>      Double_Kernel;
#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
typedef CGAL::Simple_cartesian<CGAL::Exact_algebraic>  Algebraic_Kernel;
#endif

typedef CGAL::Simple_cartesian<CGAL::Exact_rational>  Rational_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Exact_integer>  Integer_Kernel;


typedef CGAL::Integral_domain_without_division_tag        Ring;
typedef CGAL::Field_tag       Field;
typedef CGAL::Field_with_sqrt_tag  Sqrt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Double_Kernel,Ring>
Double_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Double_Kernel,Sqrt>
Double_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Double_Kernel,Field>
Double_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Double_Kernel,Sqrt>
Double_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<IT_Kernel,Ring>
IT_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<IT_Kernel,Sqrt>
IT_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<IT_Kernel,Field>
IT_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<IT_Kernel,Sqrt>
IT_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<IF_Kernel,Ring>
IF_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<IF_Kernel,Sqrt>
IF_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<IF_Kernel,Field>
IF_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<IF_Kernel,Sqrt>
IF_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Algebraic_Kernel,Ring>
Algebraic_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Algebraic_Kernel,Sqrt>
Algebraic_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Algebraic_Kernel,Field>
Algebraic_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Algebraic_Kernel,Sqrt>
Algebraic_Sqrt_Gt;
#endif

//----------------------------------------------------------------------


typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Integer_Kernel,Ring>
Integer_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Rational_Kernel,Ring>
Rational_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Rational_Kernel,Field>
Rational_Field_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Ring>
F_Algebraic_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Sqrt>
F_Algebraic_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Field>
F_Algebraic_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Sqrt>
F_Algebraic_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Integer_Kernel,Ring>
F_Integer_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Rational_Kernel,Ring>
F_Rational_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2
<Double_Kernel,Sqrt,Rational_Kernel,Field>
F_Rational_Field_Gt;


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------


template<class Gt>
typename Gt::Site_2 to_site(typename Gt::Point_2 p)
{
  return Gt::Site_2::construct_site_2(p);
}

template<class Gt>
typename Gt::Site_2 to_site(typename Gt::Segment_2 s)
{
  return Gt::Site_2::construct_site_2(s.source(), s.target());
}

template<class Gt, typename A, typename B, typename C, typename D>
void test_incircle(const A & p, const B & q, const C & r, const D & t,
    const CGAL::Sign & sign)
{
  typedef typename Gt::Vertex_conflict_2   Vertex_conflict_2;
  typedef typename Gt::Site_2   Site_2;
  Gt gt;
  Vertex_conflict_2 incircle = gt.vertex_conflict_2_object();
  Site_2 sp = to_site<Gt>(p);
  Site_2 sq = to_site<Gt>(q);
  Site_2 sr = to_site<Gt>(r);
  Site_2 st = to_site<Gt>(t);
  CGAL::Sign s = incircle(sp, sq, sr, st);
  std::cout << "test: " << sp << " " << sq << " " << sr << "  " << st;
  std::cout << "   " << sign << " " << s;
  std::cout << std::endl;
  assert(s == sign);
}

// ((p, q, inf), t) incircle test
template<class Gt, typename A, typename B, typename D>
void test_incircle(const A & p, const B & q, const D & t,
    const CGAL::Sign & sign)
{
  typedef typename Gt::Vertex_conflict_2   Vertex_conflict_2;
  typedef typename Gt::Site_2   Site_2;
  Gt gt;
  Vertex_conflict_2 incircle = gt.vertex_conflict_2_object();
  Site_2 sp = to_site<Gt>(p);
  Site_2 sq = to_site<Gt>(q);
  Site_2 st = to_site<Gt>(t);
  CGAL::Sign s = incircle(sp, sq, st);
  std::cout << "test: " << sp << " " << sq << " " << "inf" << "  " << st;
  std::cout << "   " << sign << " " << s;
  std::cout << std::endl;
  assert(s == sign);
}


template<class Gt>
void test_traits(const char* title)
{
  typedef typename Gt::Point_2             Point_2;
  typedef typename Gt::Segment_2           Segment_2;

  std::cout << "====================================" << std::endl;
  std::cout << title << std::endl;
  std::cout << "------------------------------------" << std::endl;

  // ((p,q,inf),t) tests (vertex at infinity)
  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100,0),
      Point_2(150, 30),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100,0),
      Point_2(50, -20),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100,0),
      Point_2(150, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100,0),
      Point_2(-10, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100,0),
      Point_2(40, 0),
      CGAL::NEGATIVE);


  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(20, 10),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(60, 10),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(60, 15),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(40, -30),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(-20, 10),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(60, 60),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(10, 10),
      Point_2(60, 30),
      Point_2(40, 60),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Point_2(20, 70),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Point_2(1, 69),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Point_2(70, 70),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Point_2(-10, 45),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Point_2(-5, 85),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Point_2(30, 70),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Segment_2(Point_2(20, 70), Point_2(10, 75)),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 68),
      Segment_2(Point_2(20, 70), Point_2(1, 50)),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(26, 70),
      Point_2(1, 40),
      Segment_2(Point_2(20, 70), Point_2(1, 55)),
      CGAL::NEGATIVE);

  // PPPP tests
  // PPPP three points at corners
  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(150, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(100, 100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(51, 35),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(0, 68),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(50, 0),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(79, 100),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 100),
      Point_2(0, 0),
      Point_2(100, 0),
      Point_2(100, 32),
      CGAL::NEGATIVE);

  // PPPP one point at NW corner
  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(85, -32),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(-30, 13),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(-25, -50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(75, -50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(75, 50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(50, -50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 50),
      Point_2(0, -50),
      Point_2(75, 0),
      Point_2(-25, -25),
      CGAL::NEGATIVE);

  // PPPP two points at corners (of same side)
  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(-25, -25),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(0, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(0, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(0, 75),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(0, 50),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(25, 100),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(50, 0),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(75, 70),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(100, 0),
      Point_2(100, 100),
      Point_2(0, 25),
      Point_2(100, 33),
      CGAL::NEGATIVE);

  // PPPP same side points and other point between
  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(45, 73),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(50, 50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(50, -50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(-50, 50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(-50, -50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(-10, -50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(-1, -2),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(10, 50),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(+25, 50),
      Point_2(-25, 50),
      Point_2(10, -50),
      Point_2(0, -50),
      CGAL::NEGATIVE);

  // PPPP same side and other same coordinate
  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(-5, -10),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(45, -78),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(70, 27),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(105, 27),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(0, 25),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(100, 25),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(100, -75),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(0, -75),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(0, -50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(50, -75),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(50, 25),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(3, 25),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(0, -25),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(100, -20),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(0, 0),
      Point_2(100, -50),
      Point_2(100, 0),
      Point_2(98, -70),
      CGAL::NEGATIVE);

  // PPPP same side and outside
  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(50, -2),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(50, 102),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(0, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(100, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(1, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(98, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(0, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(100, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(97, 100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(0, 99),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(100, 40),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(10, 100),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(75, 100),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(3, 100),
      Point_2(25, 0),
      Point_2(75, 0),
      Point_2(50, 0),
      CGAL::NEGATIVE);

  // PPPP same side and other at corner
  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(0, 50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(0, -50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(0, -102),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(-100, -50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(-34, -50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(-100, -25),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(-100, +25),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(-96, 2),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-100, 50),
      Point_2(0, -25),
      Point_2(0, +25),
      Point_2(-50, 50),
      CGAL::NEGATIVE);

  // PPPP same side and other close
  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(100, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(0, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(0, -100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(-2, -49),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(80, 0),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(100, -20),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(100, -70),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(100, -100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(100, -50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(50, -100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(20, -100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(48, 0),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(0, -40),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(20, -25),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(20, 0),
      Point_2(0, -25),
      Point_2(0, -75),
      Point_2(5, -75),
      CGAL::NEGATIVE);

  // PPPP same coordinate, opposite side
  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(3, 34),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-23, -2),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(0, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-100, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-100, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(0, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(0, 25),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-75, 100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-100, 75),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-100, 50),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-100, 25),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-75, 0),
      CGAL::ZERO);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(0, 50),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-50, 100),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-50, 0),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-2, 32),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(-25, 100),
      Point_2(-25, 0),
      Point_2(0, 75),
      Point_2(-99, 97),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(200, 0),
      Point_2(180, 160),
      Point_2(180, 140),
      Point_2(160, 160),
      CGAL::POSITIVE);

  test_incircle<Gt>(Point_2(400, 500),
                    Segment_2(Point_2(400, 500), Point_2(300, 500)),
                    Point_2(500, 200),
                    Segment_2(Point_2(500, 500), Point_2(400, 500)),
                    CGAL::POSITIVE);

  test_incircle<Gt>(Point_2(20, 0), Point_2(60, 0),
                    Segment_2(Point_2(0, -80), Point_2(0, 80)),
                    Point_2(40, 40),
                    CGAL::NEGATIVE);

  test_incircle<Gt>(Point_2(20, 0), Point_2(60, 0),
                    Segment_2(Point_2(0, -80), Point_2(0, 80)),
                    Point_2(40, -40),
                    CGAL::POSITIVE);

  test_incircle<Gt>(Point_2(-23, 4), Point_2(-17, 37),
                    Segment_2(Point_2(-91, 36), Point_2(36, 87)),
                    Point_2(-17, 40),
                    CGAL::POSITIVE);

  test_incircle<Gt>(Point_2(0, 100),
                    Segment_2(Point_2(0, 50), Point_2(0, 100)),
                    Segment_2(Point_2(-50, 50), Point_2(50, -50)),
                    Point_2(50, 0),
                    CGAL::NEGATIVE);

  test_incircle<Gt>(Segment_2(Point_2(-50, 50), Point_2(50, -50)),
                    Point_2(0, 100),
                    Segment_2(Point_2(0, 50), Point_2(0, 100)),
                    Point_2(50, 0),
                    CGAL::NEGATIVE);

  test_incircle<Gt>(Segment_2(Point_2(0, 50), Point_2(0, 100)),
                    Segment_2(Point_2(-50, 50), Point_2(50, -50)),
                    Point_2(0, 100),
                    Point_2(50, 0),
                    CGAL::NEGATIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(60, 40), Point_2(70, 60)),
      Point_2(100, 40),
      Segment_2(Point_2(30, 110), Point_2(100, 40)),
      Point_2(60, 20),
      CGAL::ZERO);

  test_incircle<Gt>(
      Segment_2(Point_2(60, 40), Point_2(70, 60)),
      Point_2(60, 20),
      Point_2(100, 40),
      Segment_2(Point_2(30, 110), Point_2(100, 40)),
      CGAL::ZERO);

  test_incircle<Gt>(
      Segment_2(Point_2(-100, -50), Point_2(50, 100)),
      Segment_2(Point_2(50, 50), Point_2(100, 50)),
      Point_2(50, 50),
      Point_2(0, 0),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(-100, -50), Point_2(50, 100)),
      Segment_2(Point_2(100, 0), Point_2(100, -100)),
      Point_2(100, 0),
      Point_2(0, 0),
      CGAL::ZERO);

  test_incircle<Gt>(
      Segment_2(Point_2(-100, -50), Point_2(50, 100)),
      Point_2(100, 0),
      Segment_2(Point_2(100, 0), Point_2(100, 100)),
      Point_2(0, 0),
      CGAL::ZERO);

  // 1seg1hsegnoseg.cin validity test
  test_incircle<Gt>(
      Point_2(60, 40),
      Segment_2(Point_2(10, 120), Point_2(60, 20)),
      Point_2(70, 40),
      Point_2(60, 20),
      CGAL::POSITIVE);

  // 3segstepnosegbef3.cin validity test
  test_incircle<Gt>(
      Point_2(0, 50),
      Segment_2(Point_2(-50, 50), Point_2(50, -50)),
      Point_2(0, 100),
      Point_2(-50, 50),
      CGAL::POSITIVE);

  // 2a0minimalbeforenoseg.cin validity test
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(-50, 150)),
      Point_2(0, 100),
      Point_2(100, 100),
      Point_2(0, 0),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(-100, 50), Point_2(100, 250)),
      Point_2(50, 0),
      Point_2(100, 0),
      Point_2(200, -50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(-100, 50), Point_2(100, 250)),
      Point_2(50, 0),
      Point_2(100, 0),
      Point_2(150, 100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Segment_2(Point_2(-100, 50), Point_2(100, 250)),
      Point_2(50, 0),
      Point_2(100, 0),
      Point_2(50, 100),
      CGAL::NEGATIVE);

  //Start: Tests for pps case with vertical segment and two points with same y coordinate on the right of seg
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(200, 50),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(200, 10),
      CGAL::POSITIVE);
  // t is at the bottom right corner of the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(200, -50),
      CGAL::POSITIVE);
  // t is at the top right corner of the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(200, 150),
      CGAL::POSITIVE);
  // t is on bottom horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(20, -50),
      CGAL::NEGATIVE);
  // t is on bottom horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(180, -50),
      CGAL::NEGATIVE);
  // t is on top horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(20, 150),
      CGAL::POSITIVE);
  // t is on top horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(180, 150),
      CGAL::POSITIVE);
  // t is on the supporting line of the segment
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(0, -40),
      CGAL::POSITIVE);
  // t is at the corner of the square along the segment
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(0, -50),
      CGAL::POSITIVE);
  // t is inside the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(50, 40),
      CGAL::NEGATIVE);
  // t is outside the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(200, 20),
      Point_2(200, 60),
      Point_2(20, -100),
      CGAL::POSITIVE);

  //End: Tests for pps case with vertical segment and two points with same y coordinate on the right of seg

  //Start: Tests for pps case with vertical segment and two points with same y coordinate on the left of seg
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-200, 50),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-200, 10),
      CGAL::POSITIVE);
  // t is at the bottom left corner of the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-200, -50),
      CGAL::POSITIVE);
  // t is at the top left corner of the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-200, 150),
      CGAL::POSITIVE);
  // t is on bottom horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-20, -50),
      CGAL::NEGATIVE);
  // t is on bottom horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-180, -50),
      CGAL::NEGATIVE);
  // t is on top horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-20, 150),
      CGAL::POSITIVE);
  // t is on top horizontal edge of square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-180, 150),
      CGAL::POSITIVE);
  // t is on the supporting line of the segment
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(0, -40),
      CGAL::POSITIVE);
  // t is at the corner of the square along the segment
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(0, -50),
      CGAL::POSITIVE);
  //t inside the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-50, 40),
      CGAL::NEGATIVE);
  //t outside the square
  test_incircle<Gt>(
      Segment_2(Point_2(0, 0), Point_2(0, 100)),
      Point_2(-200, 60),
      Point_2(-200, 20),
      Point_2(-20, -100),
      CGAL::POSITIVE);

  //End: Tests for pps case with vertical segment and two points with same y coordinate on the left of seg

  // PSS case vertex computation, pssprob1.cin
  test_incircle<Gt>(
      Segment_2(Point_2(0, 100), Point_2(100, 100)),
      Segment_2(Point_2(0, 100), Point_2(  0,   0)),
      Point_2(20, -40),
      Point_2(150, 30),
      CGAL::POSITIVE);

  // SSS axis-parallel segments
  test_incircle<Gt>(
      Segment_2(Point_2(200,  50), Point_2(200, 150)),
      Segment_2(Point_2(-50, 200), Point_2(450, 200)),
      Segment_2(Point_2(-50,   0), Point_2(450,   0)),
      Point_2(250, 50),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(200,  50), Point_2(200, 150)),
      Segment_2(Point_2(-50, 200), Point_2(450, 200)),
      Segment_2(Point_2(-50,   0), Point_2(450,   0)),
      Point_2(-50, 100),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(200,  50), Point_2(200, 150)),
      Segment_2(Point_2(-50, 200), Point_2(450, 200)),
      Segment_2(Point_2(-50,   0), Point_2(450,   0)),
      Point_2(0, 100),
      CGAL::ZERO);

  test_incircle<Gt>(
      Segment_2(Point_2(200,  50), Point_2(200, 150)),
      Segment_2(Point_2(-50, 200), Point_2(450, 200)),
      Segment_2(Point_2(-50,   0), Point_2(450,   0)),
      Point_2(100, 100),
      CGAL::NEGATIVE);

  // SSS sssoriented8 vertex:
  // s 100 -50 200 50
  // s 150 250 150 50
  // s -100 100 200 -200
  test_incircle<Gt>(
      Segment_2(Point_2(100, -50), Point_2(200, 50)),
      Segment_2(Point_2(150, 250), Point_2(150, 50)),
      Segment_2(Point_2(-100, 100), Point_2(200, -200)),
      Point_2(165, -32),
      CGAL::POSITIVE);

  test_incircle<Gt>(
      Segment_2(Point_2(100, -50), Point_2(200, 50)),
      Segment_2(Point_2(150, 250), Point_2(150, 50)),
      Segment_2(Point_2(-100, 100), Point_2(200, -200)),
      Point_2(110, 33),
      CGAL::NEGATIVE);

  // n31.cin related test
  test_incircle<Gt>(
      Segment_2(Point_2(-100, 50), Point_2(0, 0)),
      Segment_2(Point_2(0, 0), Point_2(60, 0)),
      Segment_2(Point_2(60, 0), Point_2(130, 40)),
      Point_2(60, 0),
      CGAL::POSITIVE);

  // tests related to points_inside_touching_sides, e.g. br83.cin
  test_incircle<Gt>(
      Point_2(80, 40),
      Segment_2(Point_2(50, 20), Point_2(70, 60)),
      Segment_2(Point_2(30, 110), Point_2(120, 20)),
      Point_2(60, 20),
      CGAL::NEGATIVE);

  test_incircle<Gt>(
      Point_2(80, 40),
      Point_2(60, 20),
      Segment_2(Point_2(30, 110), Point_2(120, 20)),
      Segment_2(Point_2(50, 20), Point_2(70, 60)),
      CGAL::POSITIVE);


  // PSS bdiff=2 point on opposite side pssd2btw3.cin

  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(0, 100),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(0, 50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(100, -50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(100, 250),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(250, 100),
       CGAL::POSITIVE);

  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(200, 100),
       CGAL::ZERO);

  test_incircle<Gt>(
       Segment_2(Point_2(100, 250), Point_2(0, 50)),
       Segment_2(Point_2(0, 50), Point_2(100, -50)),
       Point_2(200, 50),
       Point_2(100, 50),
       CGAL::NEGATIVE);


  // PSS bdiff=4 point on corner pssd4a.cin

  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(-50, -50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(0, 100),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(100, -150),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(100, 100),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(100, 150),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(200, -50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(50, 100),
       CGAL::POSITIVE);

  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Point_2(0, 50),
       Point_2(150, -100),
       CGAL::ZERO);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(0, 50),
       CGAL::ZERO);

  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(150, -50),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(50, -100),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(100, 150), Point_2(200, -50)),
       Segment_2(Point_2(-100, -50), Point_2(100, -150)),
       Point_2(150, -100),
       Point_2(12, -27),
       CGAL::NEGATIVE);

  // PSS bdiff=3 point on a side pssd3a.cin

  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(-100, -50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(-100, 0),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(-100, 100),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(-50, 150),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(0, 200),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(150, 50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(0, 0),
       Point_2(-50, 50),
       CGAL::POSITIVE);

  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(0, 100),
       CGAL::ZERO);

  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(0, 0),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(90, 90),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(90, -25),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-100, -50), Point_2(150, -50)),
       Segment_2(Point_2(0, 200), Point_2(150, 50)),
       Point_2(-50, 50),
       Point_2(-25, -35),
       CGAL::NEGATIVE);


  // PSS bdiff=3 point on a side variation pssd3varxa.cin

  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(-150, 50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(-50, 0),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(-50, 50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(150, -50),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(150, 200),
       CGAL::POSITIVE);
  test_incircle<Gt>(
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(100, 50),
       Point_2(0, 100),
       CGAL::POSITIVE);

  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(50, 150),
       CGAL::ZERO);

  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(100, 50),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(50, 50),
       CGAL::NEGATIVE);
  test_incircle<Gt>(
       Point_2(0, 100),
       Segment_2(Point_2(-150, 50), Point_2(150, -50)),
       Segment_2(Point_2(150, -50), Point_2(150, 200)),
       Point_2(90, 80),
       CGAL::NEGATIVE);

  // pssd5aless2.cin related validity test:
  // At the moment, we expect zero.
  // A negative value could also be possible (if points have priority
  // over segments), but this would require more changes in other
  // predicates.
  test_incircle<Gt>(
       Segment_2(Point_2(150, 0), Point_2(50, -100)),
       Point_2(100, 100),
       Point_2(-100, 100),
       Point_2(100, 0),
       CGAL::ZERO);

  std::cout << "====================================" << std::endl;
  std::cout << std::endl;
}


int main(int , char**)
{
#if 0
  test_traits<Double_Ring_Gtwi>("Double Ring WI");
  test_traits<Double_Sqrt_Gtwi>("Double Sqrt WI");

  test_traits<Double_Field_Gt>("Double Field");
  test_traits<Double_Sqrt_Gt>("Double Sqrt");

  std::cout << std::endl;
  std::cout << "************************************"
            << "************************************" << std::endl;
  std::cout << std::endl;

  test_traits<IT_Ring_Gtwi>("IT Ring WI");
  test_traits<IT_Sqrt_Gtwi>("IT Sqrt WI");

  test_traits<IT_Field_Gt>("IT Field");
  test_traits<IT_Sqrt_Gt>("IT Sqrt");

  std::cout << std::endl;
  std::cout << "************************************"
            << "************************************" << std::endl;
  std::cout << std::endl;

  test_traits<IF_Ring_Gtwi>("IF Ring WI");
  test_traits<IF_Sqrt_Gtwi>("IF Sqrt WI");

  test_traits<IF_Field_Gt>("IF Field");
  test_traits<IF_Sqrt_Gt>("IF Sqrt");

  std::cout << std::endl;
  std::cout << "************************************"
            << "************************************" << std::endl;
  std::cout << std::endl;
#endif

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
  test_traits<Algebraic_Ring_Gtwi>("CORE Ring WI");
  test_traits<Algebraic_Sqrt_Gtwi>("CORE Sqrt WI");

  test_traits<Algebraic_Field_Gt>("CORE Field");
  test_traits<Algebraic_Sqrt_Gt>("CORE Sqrt");
#endif

  test_traits<Integer_Ring_Gtwi>("Integer Ring WI");

  test_traits<Rational_Ring_Gtwi>("Rational Ring WI");
  test_traits<Rational_Field_Gt>("Rational Field");

  std::cout << std::endl;
  std::cout << "************************************"
            << "************************************" << std::endl;
  std::cout << std::endl;

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
  test_traits<F_Algebraic_Ring_Gtwi>("F Algebraic Ring WI");
  test_traits<F_Algebraic_Sqrt_Gtwi>("F Algebraic Sqrt WI");

  test_traits<F_Algebraic_Field_Gt>("F Algebraic Field");
  test_traits<F_Algebraic_Sqrt_Gt>("F Algebraic Sqrt");
#endif

  test_traits<F_Integer_Ring_Gtwi>("F Integer Ring WI");

  test_traits<F_Rational_Ring_Gtwi>("F Rational Ring WI");
  test_traits<F_Rational_Field_Gt>("F Rational Field");

  std::cout << std::endl;
  std::cout << "************************************"
            << "************************************" << std::endl;
  std::cout << std::endl;


  return 0;
}


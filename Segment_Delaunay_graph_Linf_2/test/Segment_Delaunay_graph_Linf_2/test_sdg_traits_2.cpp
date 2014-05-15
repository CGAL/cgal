#include <CGAL/basic.h>

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#endif

typedef CGAL::Simple_cartesian<double>      Double_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<true> >   IT_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >  IF_Kernel;
typedef CGAL::Simple_cartesian<double>      Double_Kernel;
#ifdef CGAL_USE_CORE
typedef CGAL::Simple_cartesian<CORE::Expr>  CORE_Kernel;
#endif
#ifdef CGAL_USE_GMP
typedef CGAL::Simple_cartesian<CGAL::Gmpq>  Gmpq_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Gmpz>  Gmpz_Kernel;
#endif

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

#ifdef CGAL_USE_CORE
typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<CORE_Kernel,Ring>
CORE_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<CORE_Kernel,Sqrt>
CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<CORE_Kernel,Field>
CORE_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<CORE_Kernel,Sqrt>
CORE_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

#ifdef CGAL_USE_GMP
typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Gmpz_Kernel,Ring>
Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2
<Gmpq_Kernel,Ring>
Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Gmpq_Kernel,Field>
Gmpq_Field_Gt;
#endif

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifdef CGAL_USE_CORE
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,CORE_Kernel,Ring>
F_CORE_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,CORE_Kernel,Sqrt>
F_CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2
<Double_Kernel,Sqrt,CORE_Kernel,Field>
F_CORE_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2
<Double_Kernel,Sqrt,CORE_Kernel,Sqrt>
F_CORE_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

#ifdef CGAL_USE_GMP
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Gmpz_Kernel,Ring>
F_Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Gmpq_Kernel,Ring>
F_Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2
<Double_Kernel,Sqrt,Gmpq_Kernel,Field>
F_Gmpq_Field_Gt;
#endif

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


template<class Gt>
void test_traits(const char* title)
{
  typedef typename Gt::Point_2             Point_2;
  typedef typename Gt::Segment_2           Segment_2;
  typedef typename Gt::Site_2              Site_2;

  std::cout << "====================================" << std::endl;
  std::cout << title << std::endl;
  std::cout << "------------------------------------" << std::endl;

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

#ifdef CGAL_USE_CORE
  test_traits<CORE_Ring_Gtwi>("CORE Ring WI");
  test_traits<CORE_Sqrt_Gtwi>("CORE Sqrt WI");

  test_traits<CORE_Field_Gt>("CORE Field");
  test_traits<CORE_Sqrt_Gt>("CORE Sqrt");
#endif

#ifdef CGAL_USE_GMP
  test_traits<Gmpz_Ring_Gtwi>("Gmpz Ring WI");

  test_traits<Gmpq_Ring_Gtwi>("Gmpq Ring WI");
  test_traits<Gmpq_Field_Gt>("Gmpq Field");
#endif

  std::cout << std::endl;
  std::cout << "************************************"
	    << "************************************" << std::endl;
  std::cout << std::endl;

#ifdef CGAL_USE_CORE
  test_traits<F_CORE_Ring_Gtwi>("F CORE Ring WI");
  test_traits<F_CORE_Sqrt_Gtwi>("F CORE Sqrt WI");

  test_traits<F_CORE_Field_Gt>("F CORE Field");
  test_traits<F_CORE_Sqrt_Gt>("F CORE Sqrt");
#endif

#ifdef CGAL_USE_GMP
  test_traits<F_Gmpz_Ring_Gtwi>("F Gmpz Ring WI");

  test_traits<F_Gmpq_Ring_Gtwi>("F Gmpq Ring WI");
  test_traits<F_Gmpq_Field_Gt>("F Gmpq Field");
#endif

  std::cout << std::endl;
  std::cout << "************************************"
	    << "************************************" << std::endl;
  std::cout << std::endl;


  return 0;
}

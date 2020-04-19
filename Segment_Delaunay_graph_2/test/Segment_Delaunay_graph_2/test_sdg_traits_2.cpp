#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
#include <CGAL/Exact_algebraic.h>
#endif

#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_integer.h>


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

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Double_Kernel,Ring>
Double_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Double_Kernel,Sqrt>
Double_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<Double_Kernel,Field>
Double_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_traits_2<Double_Kernel,Sqrt>
Double_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<IT_Kernel,Ring>
IT_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<IT_Kernel,Sqrt>
IT_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<IT_Kernel,Field>
IT_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_traits_2<IT_Kernel,Sqrt>
IT_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<IF_Kernel,Ring>
IF_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<IF_Kernel,Sqrt>
IF_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<IF_Kernel,Field>
IF_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_traits_2<IF_Kernel,Sqrt>
IF_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Algebraic_Kernel,Ring>
Algebraic_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Algebraic_Kernel,Sqrt>
Algebraic_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<Algebraic_Kernel,Field>
Algebraic_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_traits_2<Algebraic_Kernel,Sqrt>
Algebraic_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Integer_Kernel,Ring>
Integer_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Rational_Kernel,Ring>
Rational_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<Rational_Kernel,Field>
Rational_Field_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Ring>
F_Algebraic_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Sqrt>
F_Algebraic_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Field>
F_Algebraic_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2
<Double_Kernel,Sqrt,Algebraic_Kernel,Sqrt>
F_Algebraic_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Integer_Kernel,Ring>
F_Integer_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Rational_Kernel,Ring>
F_Rational_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2
<Double_Kernel,Sqrt,Rational_Kernel,Field>
F_Rational_Field_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<class Gt>
struct Get_sites
{
  template<class OutputIterator>
  OutputIterator operator()(OutputIterator it) const
  {
    typedef typename Gt::Point_2  Point_2;
    typedef typename Gt::Site_2   Site_2;
    Point_2 p1(400, 500);
    Point_2 p2(300, 500);
    Point_2 p3(500, 200);
    Point_2 p4(500, 500);

    Site_2 s1 = Site_2::construct_site_2(p1);
    Site_2 s2 = Site_2::construct_site_2(p1, p2);
    Site_2 s3 = Site_2::construct_site_2(p3);
    Site_2 s4 = Site_2::construct_site_2(p4, p1);

    *it++ = s1;
    *it++ = s2;
    *it++ = s3;
    *it++ = s4;
    return it;
  }
};

template<class Gt>
bool test_traits(const char* title)
{
  typedef typename Gt::Site_2              Site_2;
  typedef typename Gt::Vertex_conflict_2   Vertex_conflict_2;

  typedef Get_sites<Gt> Site_creator;

  Gt gt;
  Vertex_conflict_2 incircle = gt.vertex_conflict_2_object();
  Site_creator creator;

  std::vector<Site_2>  svec;
  creator(std::back_inserter(svec));

  std::cout << "====================================" << std::endl;
  std::cout << title << std::endl;
  std::cout << "------------------------------------" << std::endl
            << std::endl;
  std::cout << "Sites: " << std::endl;
  for (unsigned int i = 0; i < svec.size(); i++) {
    std::cout << "   " << svec[i] << std::endl;
  }
  CGAL::Sign s = incircle(svec[0], svec[1], svec[2], svec[3]);
  std::cout << "    incircle: " << s << std::endl;
  std::cout << "====================================" << std::endl;

  std::cout << std::endl;

  assert( s == CGAL::ZERO );
  return s == CGAL::ZERO;
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
  test_traits<Algebraic_Ring_Gtwi>("Algebraic Ring WI");
  test_traits<Algebraic_Sqrt_Gtwi>("Algebraic Sqrt WI");

  test_traits<Algebraic_Field_Gt>("Algebraic Field");
  test_traits<Algebraic_Sqrt_Gt>("Algebraic Sqrt");
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

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

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

#ifdef CGAL_USE_CORE
typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<CORE_Kernel,Ring>
CORE_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<CORE_Kernel,Sqrt>
CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<CORE_Kernel,Field>
CORE_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_traits_2<CORE_Kernel,Sqrt>
CORE_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

#ifdef CGAL_USE_GMP
typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Gmpz_Kernel,Ring>
Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_traits_without_intersections_2
<Gmpq_Kernel,Ring>
Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_traits_2<Gmpq_Kernel,Field>
Gmpq_Field_Gt;
#endif

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifdef CGAL_USE_CORE
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,CORE_Kernel,Ring>
F_CORE_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,CORE_Kernel,Sqrt>
F_CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2
<Double_Kernel,Sqrt,CORE_Kernel,Field>
F_CORE_Field_Gt;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2
<Double_Kernel,Sqrt,CORE_Kernel,Sqrt>
F_CORE_Sqrt_Gt;
#endif

//----------------------------------------------------------------------

#ifdef CGAL_USE_GMP
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Gmpz_Kernel,Ring>
F_Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Gmpq_Kernel,Ring>
F_Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2
<Double_Kernel,Sqrt,Gmpq_Kernel,Field>
F_Gmpq_Field_Gt;
#endif

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

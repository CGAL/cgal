#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

#include <CGAL/Filtered_exact.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>
#include <CGAL/Number_type_traits.h>

#include <CGAL/CORE_Expr.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>

typedef CGAL::Simple_cartesian<double>      Double_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<true> >   IT_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >  IF_Kernel;
typedef CGAL::Simple_cartesian<double>      Double_Kernel;
typedef CGAL::Simple_cartesian<CORE::Expr>  CORE_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Gmpq>  Gmpq_Kernel;
typedef CGAL::Simple_cartesian<CGAL::Gmpz>  Gmpz_Kernel;

template<class ET>
struct FE_Kernel
  : public CGAL::Simple_cartesian< CGAL::Filtered_exact<double,ET> >
{};

typedef FE_Kernel<CORE::Expr> FE_CORE_Kernel;
typedef FE_Kernel<CGAL::Gmpz> FE_Gmpz_Kernel;
typedef FE_Kernel<CGAL::Gmpq> FE_Gmpq_Kernel;

typedef CGAL::Ring_tag        Ring;
typedef CGAL::Field_tag       Field;
typedef CGAL::Sqrt_field_tag  Sqrt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<Double_Kernel,Ring>
Double_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<Double_Kernel,Sqrt>
Double_Sqrt_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<Double_Kernel,Field>
Double_Field_Gt;

typedef CGAL::Segment_Voronoi_diagram_traits_2<Double_Kernel,Sqrt>
Double_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<IT_Kernel,Ring>
IT_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<IT_Kernel,Sqrt>
IT_Sqrt_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<IT_Kernel,Field>
IT_Field_Gt;

typedef CGAL::Segment_Voronoi_diagram_traits_2<IT_Kernel,Sqrt>
IT_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<IF_Kernel,Ring>
IF_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<IF_Kernel,Sqrt>
IF_Sqrt_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<IF_Kernel,Field>
IF_Field_Gt;

typedef CGAL::Segment_Voronoi_diagram_traits_2<IF_Kernel,Sqrt>
IF_Sqrt_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<CORE_Kernel,Ring>
CORE_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<CORE_Kernel,Sqrt>
CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<CORE_Kernel,Field>
CORE_Field_Gt;

typedef CGAL::Segment_Voronoi_diagram_traits_2<CORE_Kernel,Sqrt>
CORE_Sqrt_Gt;

//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<Gmpz_Kernel,Ring>
Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<Gmpq_Kernel,Ring>
Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<Gmpq_Kernel,Field>
Gmpq_Field_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,CORE_Kernel,Ring>
F_CORE_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,CORE_Kernel,Sqrt>
F_CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_2
<Double_Kernel,Sqrt,CORE_Kernel,Field>
F_CORE_Field_Gt;

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_2
<Double_Kernel,Sqrt,CORE_Kernel,Sqrt>
F_CORE_Sqrt_Gt;

//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Gmpz_Kernel,Ring>
F_Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2
<Double_Kernel,Sqrt,Gmpq_Kernel,Ring>
F_Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_filtered_traits_2
<Double_Kernel,Sqrt,Gmpq_Kernel,Field>
F_Gmpq_Field_Gt;

//----------------------------------------------------------------------
//----------------------------------------------------------------------


typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<FE_CORE_Kernel,Ring>
FE_CORE_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<FE_CORE_Kernel,Sqrt>
FE_CORE_Sqrt_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<FE_CORE_Kernel,Field>
FE_CORE_Field_Gt;

typedef CGAL::Segment_Voronoi_diagram_traits_2<FE_CORE_Kernel,Sqrt>
FE_CORE_Sqrt_Gt;

//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<FE_Gmpz_Kernel,Ring>
FE_Gmpz_Ring_Gtwi;

//----------------------------------------------------------------------

typedef CGAL::Segment_Voronoi_diagram_traits_without_intersections_2
<FE_Gmpq_Kernel,Ring>
FE_Gmpq_Ring_Gtwi;

typedef CGAL::Segment_Voronoi_diagram_traits_2<FE_Gmpq_Kernel,Field>
FE_Gmpq_Field_Gt;


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace CGAL {
  unsigned int get_failures() {
    return num_failures_vertex_conflict;
  }
} // namespace CGAL

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
bool test_traits(char* title)
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
  std::cout << "------------------------------------" << std::endl;
  std::cout << "# of filter failures: " << CGAL::get_failures()
	    << std::endl << std::endl;
  std::cout << "====================================" << std::endl;
  
  std::cout << std::endl;  

  CGAL_assertion( s == CGAL::ZERO );
  return s == CGAL::ZERO;
}


int main(int argc, char* argv[])
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
  test_traits<CORE_Ring_Gtwi>("CORE Ring WI");
  test_traits<CORE_Sqrt_Gtwi>("CORE Sqrt WI");

  test_traits<CORE_Field_Gt>("CORE Field");
  test_traits<CORE_Sqrt_Gt>("CORE Sqrt");

  test_traits<Gmpz_Ring_Gtwi>("Gmpz Ring WI");

  test_traits<Gmpq_Ring_Gtwi>("Gmpq Ring WI");
  test_traits<Gmpq_Field_Gt>("Gmpq Field");

  std::cout << std::endl;
  std::cout << "************************************"
	    << "************************************" << std::endl;
  std::cout << std::endl;

  test_traits<F_CORE_Ring_Gtwi>("F CORE Ring WI");
  test_traits<F_CORE_Sqrt_Gtwi>("F CORE Sqrt WI");

  test_traits<F_CORE_Field_Gt>("F CORE Field");
  test_traits<F_CORE_Sqrt_Gt>("F CORE Sqrt");

  test_traits<F_Gmpz_Ring_Gtwi>("F Gmpz Ring WI");

  test_traits<F_Gmpq_Ring_Gtwi>("F Gmpq Ring WI");
  test_traits<F_Gmpq_Field_Gt>("F Gmpq Field");

  std::cout << std::endl;
  std::cout << "************************************"
	    << "************************************" << std::endl;
  std::cout << std::endl;

  test_traits<FE_CORE_Ring_Gtwi>("FE CORE Ring WI");
  test_traits<FE_CORE_Sqrt_Gtwi>("FE CORE Sqrt WI");

  test_traits<FE_CORE_Field_Gt>("FE CORE Field");
  test_traits<FE_CORE_Sqrt_Gt>("FE CORE Sqrt");

  test_traits<FE_Gmpz_Ring_Gtwi>("FE Gmpz Ring WI");

  test_traits<FE_Gmpq_Ring_Gtwi>("FE Gmpq Ring WI");
  test_traits<FE_Gmpq_Field_Gt>("FE Gmpq Field");

  return 0;
}

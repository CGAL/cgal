//#define CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES

#include <CGAL/Cartesian.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Lazy_circular_kernel_2.h>
#include <CGAL/Filtered_bbox_circular_kernel_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>
#include <CGAL/Timer.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <CGAL/IO/Dxf_variant_reader.h>
#include <fstream>
#include <iomanip>

#include <boost/type_traits.hpp>

// CIRCULAR KERNEL TYPEDEFS
//typedef CGAL::MP_Float RT;
//typedef CGAL::Quotient<RT> NT1;
typedef CGAL::Gmpz RT;
typedef CGAL::Gmpq NT1;
typedef CGAL::Cartesian<NT1> Linear_k1;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1> Algebraic_k1;
typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1> CircularKernel;
typedef CGAL::Arr_circular_arc_traits_2<CircularKernel> CircularK_CA_Traits;
typedef CircularKernel::Circular_arc_2 CircularKArc;
typedef std::vector<CircularKArc> CircularKArcContainer;
typedef CircularKernel::Circular_arc_2 Circular_arc_2;
typedef CircularKernel::Line_arc_2 Line_arc_2;
typedef CGAL::Arr_circular_line_arc_traits_2<CircularKernel>   CircularK_Variant_Traits;
typedef boost::variant< Circular_arc_2, Line_arc_2 > CircularKVarArc;
typedef std::vector<CircularKVarArc> CircularKVarArcContainer;

// LAZY KERNEL TYPEDEFS
typedef CGAL::Interval_nt_advanced NT3;
typedef CGAL::Cartesian<NT3> Linear_k3;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT3> Algebraic_k3;
typedef CGAL::Circular_kernel_2 <Linear_k3,Algebraic_k3> CK3_;
typedef CGAL::Lazy_circular_kernel_2<CircularKernel,CK3_> LazyCurvedK;
typedef CGAL::Arr_circular_arc_traits_2<LazyCurvedK> LazyCurvedK_CA_Traits;
typedef LazyCurvedK::Circular_arc_2 LazyArc;
typedef std::vector<LazyArc> LazyArcContainer;
typedef LazyCurvedK::Circular_arc_2 Circular_arc_3;
typedef LazyCurvedK::Line_arc_2 Line_arc_3;
typedef boost::variant<Circular_arc_3,Line_arc_3 > LazyVarArc;
typedef std::vector<LazyVarArc> LazyVarContainer;
//~ typedef CGAL::Arr_circular_line_arc_traits_2<LazyCurvedK,Line_arc_3,Circular_arc_3> LazyCurvedK_Variant_Traits;
typedef CGAL::Arr_circular_line_arc_traits_2<LazyCurvedK> LazyCurvedK_Variant_Traits;

// BBOX TYPEDEFS
typedef CGAL::Filtered_bbox_circular_kernel_2<CircularKernel>
  BBCircularKernel ;
typedef CGAL::Arr_circular_arc_traits_2<BBCircularKernel>
  BBCircularKernel_CA_Traits;
typedef BBCircularKernel::Circular_arc_2
  BBCircularKernelArc;
typedef std::vector<BBCircularKernelArc>
  BBCircularKernelArcContainer;
typedef BBCircularKernel::Circular_arc_2
  Circular_arc_6;
typedef BBCircularKernel::Line_arc_2
  Line_arc_6;
typedef boost::variant<Circular_arc_6,Line_arc_6 >
  BBCircVarArc;
typedef std::vector<BBCircVarArc>
  BBCircVarContainer;
typedef CGAL::Arr_circular_line_arc_traits_2<BBCircularKernel>  BBCircVariantTraits;

// BBOX(LAZY)
typedef CGAL::Filtered_bbox_circular_kernel_2<LazyCurvedK>
  BBLazyKernel ;
typedef CGAL::Arr_circular_arc_traits_2<BBLazyKernel>
  BBLazyKernel_CA_Traits;
typedef BBLazyKernel::Circular_arc_2
  BBLazyKernelArc;
typedef std::vector<BBLazyKernelArc>
  BBLazyKernelArcContainer;
typedef BBLazyKernel::Circular_arc_2
  Circular_arc_lazybb;
typedef BBLazyKernel::Line_arc_2
  Line_arc_lazybb;
typedef boost::variant<Circular_arc_lazybb,Line_arc_lazybb >
  BBLazyVarArc;
typedef std::vector<BBLazyVarArc>
  BBLazyVarContainer;
typedef CGAL::Arr_circular_line_arc_traits_2<BBLazyKernel>  BBLazyVariantTraits;

template <class CK,class Traits,class ArcContainer>
void do_main(const char *s) {

  // TYPEDEFS
  typedef typename CK::Circular_arc_2      C2;
  typedef typename CK::Line_arc_2          L2;
  typedef typename CGAL::Arrangement_2<Traits>          Pmwx;
  typedef typename CGAL::Arr_naive_point_location<Pmwx> Point_location;

  // LOADING CURVES
  ArcContainer ac;
  std::ifstream fin;
  fin.open (s);
  CGAL::variant_load<CK, C2, L2>(
    fin, std::back_inserter(ac));
  fin.close();

  std::cout << "Size:" << ac.size() << std::endl;

  // BENCHMARKING
  Pmwx _pm;
  Point_location _pl(_pm);
  struct rusage before, after;
  struct timeval utime, stime;
  getrusage(RUSAGE_SELF,&before);
  insert(_pm,ac.begin(),ac.end(),boost::false_type());
  getrusage(RUSAGE_SELF,&after);
  timersub(&(after.ru_utime),&(before.ru_utime),&utime);
  timersub(&(after.ru_stime),&(before.ru_stime),&stime);
  std::cout<<"Time="<< utime.tv_sec<<"."<< std::setw(6) <<
  std::setfill('0')<< utime.tv_usec <<std::endl;

  std::cerr << utime.tv_sec << "." << std::setw(6) <<
  std::setfill('0')<< utime.tv_usec << std::endl;

  std::cout << "The arrangement size:" << std::endl
            << "   V = " << _pm.number_of_vertices()
            << ",  E = " << _pm.number_of_edges()
            << ",  F = " << _pm.number_of_faces() << std::endl;

}

template <class CK,class Traits,class ArcContainer>
void do_main(int k) {

  // TYPEDEFS
  typedef typename CK::Circular_arc_2      C2;
  typedef typename CK::Point_2 Point_2;
  typedef typename CK::FT FT;
  typedef typename CGAL::Arrangement_2<Traits>          Pmwx;
  typedef typename CGAL::Arr_naive_point_location<Pmwx> Point_location;

  // LOADING CURVES
  ArcContainer ac;

  double cx, cy;
  const FT rft(1.0);

  // DENSE
  if(k == 0) {
    for(cx = 0.0; cx <= 10.0; cx += 0.5) {
      for(cy = 0.0; cy <= 10.0; cy += 0.5) {
        ac.push_back(typename CK::Circular_arc_2 ( typename CK::Circle_2(Point_2(cx,cy),rft) ) );
      }
    }
  }

  // VERY DENSE
  if(k == 1) {
    for(cx = 0.0; cx <= 0.2; cx += 0.01) {
      for(cy = 0.0; cy <= 0.2; cy += 0.01) {
        ac.push_back(typename CK::Circular_arc_2 ( typename CK::Circle_2(Point_2(cx,cy),rft) ));
      }
    }
  }

  // ONE CIRCLE
  if(k == 2) {
    ac.push_back(typename CK::Circular_arc_2 ( typename CK::Circle_2(Point_2(0,0),5) ));
  }

  // RANDOM CASE
  if(k == 3) {
    CGAL::Random generatorOfgenerator;
    int random_seed = generatorOfgenerator.get_int(0, 123456);
    CGAL::Random theRandom(random_seed);
    for(int i=0; i<100; i++) {
      double x = theRandom.get_double(0.0,1.0);
      double y = theRandom.get_double(0.0,1.0);
      double r = theRandom.get_double(0.00001,1.0);
      ac.push_back(typename CK::Circular_arc_2 ( typename CK::Circle_2(Point_2(x,y),FT(r)) ));
    }
  }

  // SPARSE
  if(k == 4) {
    double h = (std::sqrt(3.0)+0.01);
    for(cx = 0.0; cx <= 40.4; cx += 2.01) {
      for(cy = h; cy <= h*20; cy += h*2) {
        ac.push_back(typename CK::Circular_arc_2 ( typename CK::Circle_2(Point_2(cx+1.0,cy),rft) ));
      }
      for(cy = 0.0; cy <= h*20.0; cy += h*2) {
        ac.push_back(typename CK::Circular_arc_2 ( typename CK::Circle_2(Point_2(cx,cy),rft)));
      }
    }
  }

  std::cout << "Size:" << ac.size() << std::endl;

  // BENCHMARKING
  Pmwx _pm;
  Point_location _pl(_pm);
  struct rusage before, after;
  struct timeval utime, stime;
  getrusage(RUSAGE_SELF,&before);
  insert(_pm,ac.begin(),ac.end(),boost::false_type());
  getrusage(RUSAGE_SELF,&after);
  timersub(&(after.ru_utime),&(before.ru_utime),&utime);
  timersub(&(after.ru_stime),&(before.ru_stime),&stime);
  std::cout<<"Time="<< utime.tv_sec<<"."<< std::setw(6) <<
  std::setfill('0')<< utime.tv_usec << std::endl;

  std::cerr << utime.tv_sec << "." << std::setw(6) <<
  std::setfill('0')<< utime.tv_usec << std::endl;

  std::cout << "The arrangement size:" << std::endl
            << "   V = " << _pm.number_of_vertices()
            << ",  E = " << _pm.number_of_edges()
            << ",  F = " << _pm.number_of_faces() << std::endl;
}

int main(int argc, char* argv[]){

  const char* dxf_filename[] = { "DXF/51.dxf",
                                 "DXF/cad_l1.dxf",
                                 "DXF/cad_l2.dxf",
                                 "DXF/che_mod1.dxf",
                                 "DXF/CIOnZDraw.dxf",
                                 "DXF/mask1.dxf",
                                 "DXF/elekonta.dxf",
                                 "DXF/netlist_signal_1.dxf",
                                 "DXF/painttrack.dxf" };
  if(argc == 3) {
    int i = argv[1][0]-'0';
    int j = argv[2][0]-'0';
    if((j >= 0 && j < 9)) {
      if(i == 1) do_main<BBCircularKernel,BBCircVariantTraits, BBCircVarContainer>(dxf_filename[j]);
      if(i == 2) do_main<LazyCurvedK,LazyCurvedK_Variant_Traits, LazyVarContainer>(dxf_filename[j]);
      if(i == 3) do_main<CircularKernel,CircularK_Variant_Traits, CircularKVarArcContainer>(dxf_filename[j]);
      if(i == 4) do_main<BBLazyKernel,BBLazyVariantTraits, BBLazyVarContainer>(dxf_filename[j]);
      if((i >= 5) || (i <= 0)) std::cout << "INVALID PARAMETERS" << std::endl;
    } else {
      int k = -1;
      if(j == 9) k = 0;
      if(j == ('a'-'0')) k = 1;
      if(j == ('b'-'0')) k = 2;
      if(j == ('c'-'0')) k = 3;
      if(j == ('d'-'0')) k = 4;
      if(i == 1) do_main<BBCircularKernel,BBCircVariantTraits, BBCircVarContainer>(k);
      if(i == 2) do_main<LazyCurvedK,LazyCurvedK_Variant_Traits, LazyVarContainer>(k);
      if(i == 3) do_main<CircularKernel,CircularK_Variant_Traits, CircularKVarArcContainer>(k);
      if(i == 4) do_main<BBLazyKernel,BBLazyVariantTraits, BBLazyVarContainer>(k);
      if(i == 5) do_main<BBCircularKernel,BBCircularKernel_CA_Traits, BBCircularKernelArcContainer>(k);
      if(i == 6) do_main<LazyCurvedK,LazyCurvedK_CA_Traits, LazyArcContainer>(k);
      if(i == 7) do_main<CircularKernel,CircularK_CA_Traits, CircularKArcContainer>(k);
      if(i == 8) do_main<BBLazyKernel,BBLazyKernel_CA_Traits, BBLazyKernelArcContainer>(k);
    }
  } else std::cout << "INVALID PARAMETERS" << std::endl;

  return 0;
}

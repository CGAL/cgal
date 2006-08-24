//#define CGAL_MP_FLOAT_ALLOW_INEXACT

#include <CGAL/Cartesian.h>
#include <CGAL/Handle_for.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
//#include <CGAL/Gmpz.h>
//#include <CGAL/Gmpq.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Arr_circular_arc_traits.h>
#include <CGAL/Lazy_curved_kernel.h>
#include <CGAL/Filtered_hexagon_curved_kernel.h> 
#include <CGAL/Filtered_bbox_curved_kernel.h>
#include <CGAL/Filtered_interval_circular_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_circular_line_arc_traits.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Dxf_variant_reader.h>
#include <fstream>

#include <CGAL/Curved_kernel/Circular_arc_2.h>

// CURVED KERNEL TYPEDEFS
//typedef CGAL::Gmpz RT;
//typedef CGAL::Gmpq NT1;
typedef CGAL::MP_Float RT;
typedef CGAL::Quotient<RT> NT1;
typedef CGAL::Cartesian<NT1> Linear_k1;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1> Algebraic_k1;
typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1> CircularKernel;
typedef CGAL::Arr_circular_arc_traits<CircularKernel> CircularK_CA_Traits;
typedef CircularKernel::Circular_arc_2 CircularKArc;
typedef std::vector<CircularKArc> CircularKArcContainer;
typedef CircularKernel::Circular_arc_2 Circular_arc_2;
typedef CircularKernel::Line_arc_2 Line_arc_2;
typedef CGAL::Arr_circular_line_arc_traits<CircularKernel,
                  Line_arc_2,Circular_arc_2>   CircularK_Variant_Traits;
typedef boost::variant< Circular_arc_2, Line_arc_2 > CircularKVarArc;
typedef std::vector<CircularKVarArc> CircularKVarArcContainer; 

// LAZY KERNEL TYPEDEFS
typedef CGAL::Interval_nt_advanced NT3;
typedef CGAL::Cartesian<NT3> Linear_k3;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT3> Algebraic_k3;
typedef CGAL::Circular_kernel_2 <Linear_k3,Algebraic_k3> CK3_;
typedef CGAL::Lazy_curved_kernel<CircularKernel,CK3_> LazyCurvedK;
typedef CGAL::Arr_circular_arc_traits<LazyCurvedK> LazyCurvedK_CA_Traits;
typedef LazyCurvedK::Circular_arc_2 LazyArc;
typedef std::vector<LazyArc> LazyArcContainer;
typedef LazyCurvedK::Circular_arc_2 Circular_arc_3;
typedef LazyCurvedK::Line_arc_2 Line_arc_3; 
typedef boost::variant<Circular_arc_3,Line_arc_3 > LazyVarArc;
typedef std::vector<LazyVarArc> LazyVarContainer;
typedef CGAL::Arr_circular_line_arc_traits<LazyCurvedK,
         Line_arc_3,Circular_arc_3> LazyCurvedK_Variant_Traits;

// AH TYPEDEFS
typedef CGAL::Filtered_interval_circular_kernel<CircularKernel>
  AHCircularKernel;
typedef CGAL::Arr_circular_arc_traits<AHCircularKernel>
  AHCircularKernel_CA_Traits;
typedef AHCircularKernel::Circular_arc_2
  AHCircularKernelArc;
typedef std::vector<AHCircularKernelArc>
  AHCircularKernelArcContainer;
typedef AHCircularKernel::Circular_arc_2
  Circular_arc_8;
typedef AHCircularKernel::Line_arc_2
  Line_arc_8;
typedef boost::variant<Circular_arc_8,Line_arc_8 >
  AHCircVarArc;
typedef std::vector<AHCircVarArc> AHCircVarContainer;
typedef CGAL::Arr_circular_line_arc_traits<AHCircularKernel,
  Line_arc_8,Circular_arc_8>  AHCircVariantTraits; 

// BBOX TYPEDEFS
typedef CGAL::Filtered_bbox_curved_kernel<CircularKernel>                               
  BBCircularKernel ;
typedef CGAL::Arr_circular_arc_traits<BBCircularKernel>
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
typedef CGAL::Arr_circular_line_arc_traits<BBCircularKernel,
  Line_arc_6,Circular_arc_6>  BBCircVariantTraits; 

// BBOX(AH) TYPEDEFS
typedef CGAL::Filtered_bbox_curved_kernel<AHCircularKernel>                               
  BBAHCircularKernel ;
typedef CGAL::Arr_circular_arc_traits<BBAHCircularKernel>
  BBAHCircularKernel_CA_Traits;
typedef BBAHCircularKernel::Circular_arc_2
  BBAHCircularKernelArc;
typedef std::vector<BBAHCircularKernelArc>
  BBAHCircularKernelArcContainer;
typedef BBAHCircularKernel::Circular_arc_2 
  Circular_arc_ah;
typedef BBAHCircularKernel::Line_arc_2
  Line_arc_ah;
typedef boost::variant<Circular_arc_ah,Line_arc_ah >
  BBAHCircVarArc;
typedef std::vector<BBAHCircVarArc>
  BBAHCircVarContainer; 
typedef CGAL::Arr_circular_line_arc_traits<BBAHCircularKernel,
  Line_arc_ah,Circular_arc_ah>  BBAHCircVariantTraits; 

template <class CK,class Traits,class ArcContainer>
void do_main(char *s) {

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
  
  // BENCHMARKING
  Pmwx _pm;
  Point_location _pl(_pm);
  CGAL::Timer clck1;  
  double t1,t2;
  clck1.start();
  t1=clck1.time();
  insert_curves(_pm,ac.begin(),ac.end());
  t2=clck1.time();
  clck1.stop();
  std::cout<<"Time="<<(t2-t1)<<std::endl;
}

int main(){

  // Set the dxf_files correctly
  char* dxf_filename[] = { "DXF/51.dxf",
                           "DXF/cad_l1.dxf",
                           "DXF/cad_l2.dxf",
                           "DXF/che_mod1.dxf",
                           "DXF/CIOnZDraw.dxf",
                           "DXF/mask1.dxf",
                           "DXF/elekonta.dxf",
                           "DXF/netlist_signal_1.dxf",
                           "DXF/painttrack.dxf" };   	
  // Do some Benchs
  //do_main<AHCircularKernel,AHCircVariantTraits,
  //  AHCircVarContainer>(dxf_filename[7]);
  //do_main<BBCircularKernel,BBCircVariantTraits, 
  //  BBCircVarContainer>(dxf_filename[7]);
  //do_main<LazyCurvedK,LazyCurvedK_Variant_Traits, 
  //  LazyVarContainer>(dxf_filename[7]);
  //do_main<CircularKernel,CircularK_Variant_Traits, 
  //  CircularKVarArcContainer>(dxf_filename[7]);
  do_main<BBAHCircularKernel,BBAHCircVariantTraits, 
    BBAHCircVarContainer>(dxf_filename[3]);
  
  return 0;
}

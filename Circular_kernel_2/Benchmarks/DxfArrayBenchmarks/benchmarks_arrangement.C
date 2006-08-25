#define CGAL_CAST_INT

#define CIRCULAR_KERNEL_2
#define LAZY_CURVED_KERNEL_2
//#define CIRCULAR_KERNEL_2_FILTERED_HEXAGON   
//#define LAZY_CURVED_KERNEL_2_FILTERED_HEXAGON 
// #define CIRCULAR_KERNEL_2_FILTERED_BBOX
//#define LAZY_CURVED_KERNEL_2_FILTERED_BBOX

#define CGAL_CAST_INT
#include <CGAL/Cartesian.h>
#include <CGAL/Handle_for.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/MP_Float.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>

#include <CGAL/intersections.h>

#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Circular_arc_traits.h>
#include <CGAL/Circular_arc_traits_tracer.h>

#include <CGAL/Lazy_circular_kernel_2.h>

#include <CGAL/Filtered_hexagon_circular_kernel_2.h>

#include <CGAL/Filtered_bbox_circular_kernel_2.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Variant_traits.h>


#include <CGAL/Random.h>
#include <fstream>

#include "benchmark.h"

int main(int argc, char* argv[])
{   typedef std::list<char*> Dxffilenames;
	Dxffilenames dxffilenames;
    char* Dxffilename[]={"myFirst.dxf","minimask0.dxf","minimask1.dxf","mask0.dxf","mask1.dxf","smallpainttrack.dxf","mask0_25.dxf","mask0_5.dxf","cad_l2.dxf","cad_l1.dxf","CIOnZDraw.dxf","che_mod1.dxf","elekonta.dxf","painttrack.dxf","netlist_signal_1.dxf","51.dxf"}; 
    char* Htmlfilename;
    char* Texfilename;
    int i;
    i=0;

		
	Htmlfilename="benchmarks.html";
	Texfilename="benchmarks.tex";


for(int n = 0; n < 3; n++){
	dxffilenames.push_back(Dxffilename[n]);
	};


 Bench bench(dxffilenames,Htmlfilename,Texfilename);


for(i=1;i<6;i++){
 


/*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Circular_Kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
  #ifdef CIRCULAR_KERNEL_2
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>                      Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1>         CircularKernel;

  typedef CircularKernel::Circular_arc_2                                  Circular_arc_2;
  typedef CircularKernel::Line_arc_2                                      Line_arc_2;
  typedef CGAL::Variant_traits<CircularKernel,Line_arc_2,Circular_arc_2>  CircularK_Variant_Traits;
 
  typedef boost::variant< Circular_arc_2, Line_arc_2 >        CircularKVarArc;
  typedef std::vector<CircularKVarArc>                        CircularKVarArcContainer; 
  
  bench.kernel("CK VarTraits");
  
  bench.ComputeArrayDxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(dxffilenames);
  
//bench.Compute_dxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(Dxffilename[i]);
    
 #endif
/*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Lazy_curved_Kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
    #ifdef LAZY_CURVED_KERNEL_2
//  typedef CGAL::Quotient<CGAL::MP_Float>                       NT2;
  typedef CGAL::Gmpq                       NT2;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT2>                      Algebraic_k2;
  typedef CGAL::Circular_kernel_2 <Linear_k2, Algebraic_k2>         CK2_;
  

  typedef CGAL::Interval_nt_advanced                          NT3;
  typedef CGAL::Cartesian<NT3>                                 Linear_k3;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT3>                      Algebraic_k3;
  typedef CGAL::Circular_kernel_2 <Linear_k3,Algebraic_k3>          CK3_;
  

  typedef CGAL::Lazy_circular_kernel_2<CK2_,CK3_>                  LazyCurvedK;
  

 
  typedef LazyCurvedK::Circular_arc_2  Circular_arc_3;
  typedef LazyCurvedK::Line_arc_2  Line_arc_3; 
  typedef boost::variant<Circular_arc_3,Line_arc_3 >               LazyVarArc;
  typedef std::vector<LazyVarArc>                                  LazyVarContainer;
  typedef CGAL::Variant_traits<LazyCurvedK,Line_arc_3,Circular_arc_3> LazyCurvedK_Variant_Traits;
  
  bench.kernel("LazyK.  VarTraits");
  
   bench.ComputeArrayDxf<LazyCurvedK,LazyCurvedK_Variant_Traits,LazyVarContainer>(dxffilenames);
   //bench.Compute_dxf<LazyCurvedK,LazyCurvedK_Variant_Traits,LazyVarContainer>(Dxffilename[i]);
    
 #endif
  /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Filtered_hexagone_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
  #ifdef CIRCULAR_KERNEL_2_FILTERED_HEXAGON  

  typedef CGAL::Filtered_hexagon_circular_kernel_2<CircularKernel>        CircularKernelHexagon;

  typedef CircularKernelHexagon::Circular_arc_2                                                   Circular_arc_4;
  typedef CircularKernelHexagon::Line_arc_2                                                       Line_arc_4;
  typedef boost::variant<  Circular_arc_4, Line_arc_4 >                          CircularKernHexVarArc;
  typedef std::vector<CircularKernHexVarArc>                                     CircularKernHexVarArcContainer; 
  typedef CGAL::Variant_traits<CircularKernelHexagon,Circular_arc_4,Line_arc_4>  CircularKernHex_Variant_Traits;
  
  bench.kernel("CK Hex VarTraits");
 
 // bench.Compute<CircularKernelHexagon,CircularKernHex_Variant_Traits,CircularKernHexVarArcContainer>(Dxffilename[i]);
   bench.Compute_no_dxf<CircularKernelHexagon,CircularKernHex_Variant_Traits,CircularKernHexVarArcContainer>();

  
 #endif
  /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Filtered_hexagone_Lazy_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
 #ifdef LAZY_CURVED_KERNEL_2_FILTERED_HEXAGON 

  typedef CGAL::Filtered_hexagon_curved_kernel<LazyCurvedK>  LazyKernelHexagon;	
  
  typedef LazyKernelHexagon::Circular_arc_2                                        Circular_arc_5;
  typedef LazyKernelHexagon::Line_arc_2                                            Line_arc_5;
  typedef boost::variant<Circular_arc_5,Line_arc_5 >                  HxLazyVarArc;
  typedef std::vector<HxLazyVarArc>                                   HxLazyVarContainer;
  typedef CGAL::Variant_traits<LazyKernelHexagon,Line_arc_5,Circular_arc_5>  HxLazyVariantTraits; 
  
  bench.kernel("LazyK Hex  VarTraits") ;

  //bench.Compute<LazyKernelHexagon,HxLazyVariantTraits,HxLazyVarContainer>(Dxffilename[i]);
bench.Compute_no_dxf<LazyKernelHexagon,HxLazyVariantTraits,HxLazyVarContainer>();
 #endif
 /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!bbox_filtered_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/  
 #ifdef CIRCULAR_KERNEL_2_FILTERED_BBOX      

  typedef CGAL::Filtered_bbox_circular_kernel_2<CircularKernel>           BBCircularKernel ;
 
  typedef BBCircularKernel::Circular_arc_2                                        Circular_arc_6;
  typedef BBCircularKernel::Line_arc_2                                            Line_arc_6;
  typedef boost::variant<Circular_arc_6,Line_arc_6 >                  BBCircVarArc;
  typedef std::vector<BBCircVarArc>                                   BBCircVarContainer;
  typedef CGAL::Variant_traits<BBCircularKernel,Line_arc_6,Circular_arc_6>  BBCircVariantTraits; 
  
  bench.kernel("CK BBox VarTraits") ;
  
   bench.Compute<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>(Dxffilename[i]);
  //  bench.Compute_no_dxf<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>();
  
 #endif
 /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!bbox_hexagone_Lazy_Circular_kernel!!!!!!!!!!!!!!!!!!
						
 -------------------------------------------------------------------------------------------------------------------------*/
 #ifdef LAZY_CURVED_KERNEL_2_FILTERED_BBOX    

   typedef CGAL::Filtered_bbox_curved_kernel<LazyCurvedK>              BBLazyCurvedK; 

  typedef BBLazyCurvedK::Circular_arc_2                                        Circular_arc_7;
  typedef BBLazyCurvedK::Line_arc_2                                            Line_arc_7;
  typedef boost::variant<Circular_arc_7,Line_arc_7 >                  BBLazyVarArc;
  typedef std::vector< BBLazyVarArc>                                    BBLazyVarContainer;
  typedef CGAL::Variant_traits<BBLazyCurvedK,Line_arc_7,Circular_arc_7>   BBLazyVariantTraits; 
  
  bench.kernel("LLazyK BBox VarTraits") ;
  
  bench.Compute<BBLazyCurvedK, BBLazyVariantTraits, BBLazyVarContainer>(Dxffilename[i]);
   //bench.Compute_no_dxf<BBLazyCurvedK, BBLazyVariantTraits, BBLazyVarContainer>();
    
    
 #endif
    /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!bbox_filtered_Filtered_hexagone_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/  
       /* 
	typedef CGAL::Filtered_bbox_curved_kernel<CircularKernelHexagon>           BBCircKHexagon ;

   #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Circular_arc_traits<BBCircKHexagon>                  BBCircKHexagonCATraits;
  #else
  typedef CGAL::Circular_arc_traits<BBCircKHexagon>                  Traits0_7;
  typedef CGAL::Circular_arc_traits_tracer<Traits0_7>            BBCircKHexagonCATraits;
  #endif  
  typedef BBCircKHexagon::Circular_arc_2                              BBCircKHexagonArc;
  typedef std::vector<BBCircKHexagon>                                 BBCircKHexagonArcCont;
  bench.kernel("BBox Circular kernel filtered Hexagon CircArcTraits");

  bench.Compute_no_dxf<BBCircKHexagon,BBCircKHexagonCATraits, BBCircKHexagonArcCont>();

  typedef BBCircularKernelHexagon::Circular_arc_2                                        Circular_arc_8;
  typedef BBCircularKernelHexagon::Line_arc_2                                            Line_arc_8;
  typedef boost::variant<Circular_arc_8,Line_arc_8 >                 BBCircularKernelHexagonVarArc;
  typedef std::vector<BBCircularKernelHexagonVarArc>                                   BBCircularKernelHexagonVarContainer;
  typedef CGAL::Variant_traits<BBCircularKernelHexagon,Line_arc_8,Circular_arc_8>  BBCircularKernelHexagonVariantTraits; 
  
  bench.kernel("BBox Circular kernel filtered Hexagon  VarTraits") ;
  
   //bench.Compute<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>(Dxffilename[i]);
    bench.Compute_no_dxf<BBCircularKernelHexagon,BBCircularKernelHexagonVariantTraits,BBCircularKernelHexagonVarContainer>();*/
    /*--------------------------------------------------------------------------------------------------------------------------
  -----------------------------------------------------------------------------------------------------------------------------*/  
dxffilenames.erase(dxffilenames.begin(),dxffilenames.end());
if (i<5){
for(int n = 3*(i+1)-3; n < 3*(i+1); n++){
	dxffilenames.push_back(Dxffilename[n]);
	};
 bench.newDxfArray(dxffilenames);
 }
 }

  return 0;
};



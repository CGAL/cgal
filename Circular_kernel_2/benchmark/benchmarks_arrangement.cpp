#define CGAL_CAST_INT
#define CIRCULAR_KERNEL_2
#define LAZY_CURVED_KERNEL_2
#define CIRCULAR_KERNEL_2_FILTERED_HEXAGON   
#define LAZY_CURVED_KERNEL_2_FILTERED_HEXAGON 
#define CIRCULAR_KERNEL_2_FILTERED_BBOX
#define LAZY_CURVED_KERNEL_2_FILTERED_BBOX


#include <CGAL/Cartesian.h>
#include <CGAL/Handle_for.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/MP_Float.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>

#include <CGAL/intersections.h>

#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
//#include <CGAL/Circular_arc_traits_tracer.h>

#include <CGAL/Lazy_circular_kernel_2.h>

#ifdef CIRCULAR_KERNEL_2_FILTERED_HEXAGON
#include <CGAL/Filtered_hexagon_circular_kernel_2.h>
#endif

#include <CGAL/Filtered_bbox_circular_kernel_2.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>


#include <CGAL/Random.h>
#include <fstream>

#include "benchmark.h"

int main(int argc, char* argv[])
{	//char* Dxffilename[]={"myFirst.dxf","cad_l1.dxf"};
    const char* Dxffilename[]={"myFirst.dxf","minimask0.dxf","minimask1.dxf","mask0.dxf","mask1.dxf","smallpainttrack.dxf","mask0_25.dxf","mask0_5.dxf","cad_l2.dxf","cad_l1.dxf","CIOnZDraw.dxf","che_mod1.dxf","elekonta.dxf","painttrack.dxf","netlist_signal_1.dxf","51.dxf"}; 
    std::string Htmlfilename;
    std::string Texfilename;
    char exten[4];
    int i;
    i=0;

 
if (argv[1] != NULL)
 {
		int len =strlen(argv[1]);
		for (int j=0; j < 3 ; j++)
		{
		  exten[j]=argv[1][len - 3 + j];
		}
		if (strncmp(exten,"dxf",3) !=0)
		{
		  std::cout<< "File is not correct (*.dxf is needed)." << std::endl;
		 return 0;
		}
		else{
		Dxffilename[i] = argv[1]; 
		std::cout<< "File "<< Dxffilename[i] << " is correct."<<std::endl;
		}
		if (argc >2 and argv[2] != NULL)
		 {
			int len =strlen(argv[2]);
			for (int j=0; j < 4 ; j++)
			{
		  		exten[j]=argv[2][len - 4 + j];
			}
			if (strncmp(exten,"html",4) !=0)
			{
		 		std::cout<< "File "<< argv[2] << " is not correct (*.html is needed)." << std::endl;
		 		return 0;
			}
			else{
			std::cout<< "File "<< argv[2] << " is correct." <<std::endl;
			}
			if (argv[3] != NULL)
			{
				int len =strlen(argv[3]);
				for (int j=0; j < 3 ; j++)
				{
		  			exten[j]=argv[3][len - 3 + j];
				}
				if (strncmp(exten,"tex",3) !=0)
				{
		 			std::cout<< "File "<< argv[3] << " is not correct (*.tex is needed)."
					 << std::endl;
		 			return 0;
				}
				else{
					std::cout<< "File "<< argv[3] << " is correct." <<std::endl;
				}
 			} 
			else
	 		{	
				Texfilename="benchmarks.tex";
				std::cout<< "Default *.tex file is  :" << Texfilename<< std::endl;
 			}
 		} 
		else
 		{	
			Htmlfilename="benchmarks.html";
			Texfilename="benchmarks.tex";
			std::cout<< "Default *.html file is :" << Htmlfilename<< std::endl;
			std::cout<< "Default *.tex file is  :" << Texfilename<< std::endl;
 		}
		 
 } 
else
 {		
		Dxffilename[i]="myFirst.dxf";
 		//Dxffilename[i]="cad_l1.dxf";
		Htmlfilename="benchmarks.html";
		Texfilename="benchmarks.tex";
		std::cout<< "Default *.dxf file is  :" << Dxffilename[i]<< std::endl;
		std::cout<< "Default *.html file is :" << Htmlfilename<< std::endl;
		std::cout<< "Default *.tex file is  :" << Texfilename<< std::endl;
 }
 



 //Bench bench(Htmlfilename,Texfilename,Dxffilename[i],true);            // If you want to do benchmarks only with dxf files, you supose to use this defenition

Bench bench;			//If you want create table with all datasets you supose to use this.



// for(i=0;i<2;i++){
//  
//  std::ifstream fin;
//  fin.open (Dxffilename[i]);
//  if (!fin.is_open())
//   {
//     std::cout<<"file "<< Dxffilename[i] << " is not found"<<std::endl;
//     fin.close();
//     return 0;
// 
//   } 
//  fin.close();
/*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Circular_Kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
  #ifdef CIRCULAR_KERNEL_2

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1>     CircularKernel;

//  #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Arr_circular_arc_traits_2<CircularKernel> CircularK_CA_Traits;
//   #else
//    typedef  CGAL::Circular_arc_traits<CircularKernel>          Traits0;
//    typedef  CGAL::Circular_arc_traits_tracer<Traits0>          CircularK_CA_Traits;
//   #endif
  
  typedef CircularKernel::Circular_arc_2      CircularKArc;
  typedef std::vector<CircularKArc>           CircularKArcContainer;
  bench.kernel("CkCircArc");
  
  bench.Compute_no_dxf<CircularKernel,CircularK_CA_Traits,CircularKArcContainer>();
  
 
  typedef CircularKernel::Circular_arc_2                                  Circular_arc_2;
  typedef CircularKernel::Line_arc_2                                      Line_arc_2;
  typedef CGAL::Arr_circular_line_arc_traits_2<CircularKernel>  CircularK_Variant_Traits;
 
  typedef boost::variant< Circular_arc_2, Line_arc_2 >        CircularKVarArc;
  typedef std::vector<CircularKVarArc>                        CircularKVarArcContainer; 
  
  bench.kernel("CKVar");
  
  //bench.Compute<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(Dxffilename[i]);
  
//bench.Compute_dxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(Dxffilename[i]);
bench.Compute_no_dxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(); 

 #endif
/*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Lazy_curved_Kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
  #ifdef LAZY_CURVED_KERNEL_2

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT2;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT2>                      Algebraic_k2;
  typedef CGAL::Circular_kernel_2 <Linear_k2, Algebraic_k2>         CK2_;
  

  typedef CGAL::Interval_nt_advanced                          NT3;
  typedef CGAL::Cartesian<NT3>                                 Linear_k3;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT3>                      Algebraic_k3;
  typedef CGAL::Circular_kernel_2 <Linear_k3,Algebraic_k3>          CK3_;
  

  typedef CGAL::Lazy_circular_kernel_2<CK2_,CK3_>                  LazyCurvedK;
  
//   #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Arr_circular_arc_traits_2<LazyCurvedK> LazyCurvedK_CA_Traits;
//    #else
//    typedef  CGAL::Circular_arc_traits<LazyCurved_k>                 Traits0_2;
//    typedef  CGAL::Circular_arc_traits_tracer<Traits0_2>         LazyCurved_kTraits;
//    #endif
  
  typedef LazyCurvedK::Circular_arc_2                              LazyArc;
  typedef std::vector<LazyArc>                                  LazyArcContainer;
  
  bench.kernel("LazyCircArc") ;
  
  bench.Compute_no_dxf<LazyCurvedK,LazyCurvedK_CA_Traits,LazyArcContainer>();

 
  typedef LazyCurvedK::Circular_arc_2  Circular_arc_3;
  typedef LazyCurvedK::Line_arc_2  Line_arc_3; 
  typedef boost::variant<Circular_arc_3,Line_arc_3 >               LazyVarArc;
  typedef std::vector<LazyVarArc>                                  LazyVarContainer;
  typedef CGAL::Arr_circular_line_arc_traits_2<LazyCurvedK> LazyCurvedK_Variant_Traits;
  
  bench.kernel("LazyKVar");
  
   //bench.Compute<LazyCurvedK,LazyCurvedK_Variant_Traits,LazyVarContainer>(Dxffilename[i]);
   //bench.Compute_dxf<LazyCurvedK,LazyCurvedK_Variant_Traits,LazyVarContainer>(Dxffilename[i]);
   bench.Compute_no_dxf<LazyCurvedK,LazyCurvedK_Variant_Traits,LazyVarContainer>();

 #endif
  /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Filtered_hexagone_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
  #ifdef CIRCULAR_KERNEL_2_FILTERED_HEXAGON  

  typedef CGAL::Filtered_hexagon_circular_kernel_2<CircularKernel>        CircularKernelHexagon;
 
//   #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Arr_circular_arc_traits_2<CircularKernelHexagon>                 CircularKernHex_CA_Traits;
//   #else
//   typedef CGAL::Circular_arc_traits<CircularKernelHexagon>                 Traits0_3;
//   typedef CGAL::Circular_arc_traits_tracer<Traits0_3>     CircularKernHex_CA_Traits;
//   #endif  
  
  typedef CircularKernelHexagon::Circular_arc_2                              CircularKernHexArc;
  typedef std::vector<CircularKernHexArc>                                CircularKernHexArcContainer;
  bench.kernel("CK Hex CircArcTraits");
  
  bench.Compute_no_dxf<CircularKernelHexagon,CircularKernHex_CA_Traits,CircularKernHexArcContainer>();
 

  typedef CircularKernelHexagon::Circular_arc_2                                                   Circular_arc_4;
  typedef CircularKernelHexagon::Line_arc_2                                                       Line_arc_4;
  typedef boost::variant<  Circular_arc_4, Line_arc_4 >                          CircularKernHexVarArc;
  typedef std::vector<CircularKernHexVarArc>                                     CircularKernHexVarArcContainer; 
  typedef CGAL::Arr_circular_line_arc_traits_2<CircularKernelHexagon>  CircularKernHex_Variant_Traits;
  
  bench.kernel("CK Hex VarTraits");
 
 // bench.Compute<CircularKernelHexagon,CircularKernHex_Variant_Traits,CircularKernHexVarArcContainer>(Dxffilename[i]);
   bench.Compute_no_dxf<CircularKernelHexagon,CircularKernHex_Variant_Traits,CircularKernHexVarArcContainer>();

 
 #endif
  /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Filtered_hexagone_Lazy_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/
 #ifdef LAZY_CURVED_KERNEL_2_FILTERED_HEXAGON 


  typedef CGAL::Filtered_hexagon_circular_kernel_2<LazyCurvedK>  LazyKernelHexagon;	
   
//   #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Arr_circular_arc_traits_2<LazyKernelHexagon>                  LazyKernelHexagon_CA_Traits;
//   #else
//   typedef CGAL::Circular_arc_traits<LazyKernelHexagon>                  Traits0_4;
//   typedef CGAL::Circular_arc_traits_tracer<Traits0_4>            LazyKernelHexagon_CA_Traits;
//   #endif  

  typedef LazyKernelHexagon::Circular_arc_2                              LazyKernelHexagonArc;
  typedef std::vector<LazyKernelHexagonArc>                                  LazyKernelHexagonArcContainer;
  bench.kernel("LazyK Hex CircArcTraits");
 
  bench.Compute_no_dxf<LazyKernelHexagon,LazyKernelHexagon_CA_Traits, LazyKernelHexagonArcContainer>();
  
  typedef LazyKernelHexagon::Circular_arc_2                                        Circular_arc_5;
  typedef LazyKernelHexagon::Line_arc_2                                            Line_arc_5;
  typedef boost::variant<Circular_arc_5,Line_arc_5 >                  HxLazyVarArc;
  typedef std::vector<HxLazyVarArc>                                   HxLazyVarContainer;
  typedef CGAL::Arr_circular_line_arc_traits_2<LazyKernelHexagon>  HxLazyVariantTraits; 
  
  bench.kernel("LazyK Hex  VarTraits") ;

  //bench.Compute<LazyKernelHexagon,HxLazyVariantTraits,HxLazyVarContainer>(Dxffilename[i]);
bench.Compute_no_dxf<LazyKernelHexagon,HxLazyVariantTraits,HxLazyVarContainer>();

 #endif
 /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!bbox_filtered_Circular_kernel!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/  
 #ifdef CIRCULAR_KERNEL_2_FILTERED_BBOX   

  typedef CGAL::Filtered_bbox_circular_kernel_2<CircularKernel>           BBCircularKernel ;
 
//    #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Arr_circular_arc_traits_2<BBCircularKernel>                  BBCircularKernel_CA_Traits;
//   #else
//   typedef CGAL::Circular_arc_traits<BBCircularKernel>                  Traits0_5;
//   typedef CGAL::Circular_arc_traits_tracer<Traits0_5>            BBCircularKernel_CA_Traits;
//   #endif  
  typedef BBCircularKernel::Circular_arc_2                              BBCircularKernelArc;
  typedef std::vector<BBCircularKernelArc>                                 BBCircularKernelArcContainer;
  bench.kernel("CK BBox CircArcTraits");
 
  bench.Compute_no_dxf<BBCircularKernel,BBCircularKernel_CA_Traits, BBCircularKernelArcContainer>();
 
  typedef BBCircularKernel::Circular_arc_2                                        Circular_arc_6;
  typedef BBCircularKernel::Line_arc_2                                            Line_arc_6;
  typedef boost::variant<Circular_arc_6,Line_arc_6 >                  BBCircVarArc;
  typedef std::vector<BBCircVarArc>                                   BBCircVarContainer;
  typedef CGAL::Arr_circular_line_arc_traits_2<BBCircularKernel>  BBCircVariantTraits; 
  
  bench.kernel("CK BBox VarTraits") ;
  
  // bench.Compute<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>(Dxffilename[i]);
  bench.Compute_no_dxf<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>();

 #endif
 /*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!bbox_hexagone_Lazy_Circular_kernel!!!!!!!!!!!!!!!!!!
						
 -------------------------------------------------------------------------------------------------------------------------*/
 #ifdef LAZY_CURVED_KERNEL_2_FILTERED_BBOX  

   typedef CGAL::Filtered_bbox_circular_kernel_2<LazyCurvedK>              BBLazyCurvedK; 
  
//    #ifndef CGAL_CURVED_KERNEL_DEBUG
  typedef CGAL::Arr_circular_arc_traits_2<BBLazyCurvedK>                  BBLazyCurvedK_CA_Traits;
//   #else
//   typedef CGAL::Circular_arc_traits<BBLazyCurvedK>                  Traits0_6;
//   typedef CGAL::Circular_arc_traits_tracer<Traits0_6>            BBLazyCurvedK_CA_Traits;
//   #endif  

  typedef BBLazyCurvedK::Circular_arc_2                              BBLazyCurvedKArc;
  typedef std::vector<BBLazyCurvedKArc>                                 BBLazyCurvedKArcContainer;
  bench.kernel("LLazyK BBox CircArcTraits");
 
  bench.Compute_no_dxf<BBLazyCurvedK,BBLazyCurvedK_CA_Traits,BBLazyCurvedKArcContainer>();
 

  typedef BBLazyCurvedK::Circular_arc_2                                        Circular_arc_7;
  typedef BBLazyCurvedK::Line_arc_2                                            Line_arc_7;
  typedef boost::variant<Circular_arc_7,Line_arc_7 >                  BBLazyVarArc;
  typedef std::vector< BBLazyVarArc>                                    BBLazyVarContainer;
  typedef CGAL::Arr_circular_line_arc_traits_2<BBLazyCurvedK>   BBLazyVariantTraits; 
  
  bench.kernel("LLazyK BBox VarTraits") ;
  
  //bench.Compute<BBLazyCurvedK, BBLazyVariantTraits, BBLazyVarContainer>(Dxffilename[i]);
   bench.Compute_no_dxf<BBLazyCurvedK, BBLazyVariantTraits, BBLazyVarContainer>();
    
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
// if (i+1<2)
//   {
// try{
//  if (strcmp(Dxffilename[i+1],""))
// 	{
// 	try{	
// 	fin.open (Dxffilename[i]);
// 	}
// 	catch(...){
// 	std::cout<<"error"<<std::endl;
// 	}
//  	if (!fin.is_open())
//   		{
//     		std::cout<<"file "<< Dxffilename[i] << " is not found"<<std::endl;
// 		std::cout << "that's all" << std::endl;
//     		fin.close();
//     		break;
//   		} 
// 	else
// 		{
//  		bench.newDxfFilename(Dxffilename[i+1]);
//  		}
// 	fin.close();
// 	}
//  else
//  	{
// 	std::cout << "that's all" << std::endl;
// 	break;
//  	}
//  }
// catch(...){std::cout << "error" << std::endl;}
//  }
//  }

  return 0;
}


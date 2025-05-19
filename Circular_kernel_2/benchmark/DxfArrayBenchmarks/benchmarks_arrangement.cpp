#define CGAL_CAST_INT

#define CIRCULAR_KERNEL_2
// #define CIRCULAR_KERNEL_2_FILTERED_BBOX

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

  typedef std::variant< Circular_arc_2, Line_arc_2 >        CircularKVarArc;
  typedef std::vector<CircularKVarArc>                        CircularKVarArcContainer;

  bench.kernel("CK VarTraits");

  bench.ComputeArrayDxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(dxffilenames);

//bench.Compute_dxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(Dxffilename[i]);

 #endif

 /*-------------------------------------------------------------------------------------------------------------------------
                                                  !!!!!!!!!!!bbox_filtered_Circular_kernel!!!!!!!!!!!!!!!!!!

  -------------------------------------------------------------------------------------------------------------------------*/
 #ifdef CIRCULAR_KERNEL_2_FILTERED_BBOX

  typedef CGAL::Filtered_bbox_circular_kernel_2<CircularKernel>           BBCircularKernel ;

  typedef BBCircularKernel::Circular_arc_2                                        Circular_arc_6;
  typedef BBCircularKernel::Line_arc_2                                            Line_arc_6;
  typedef std::variant<Circular_arc_6,Line_arc_6 >                  BBCircVarArc;
  typedef std::vector<BBCircVarArc>                                   BBCircVarContainer;
  typedef CGAL::Variant_traits<BBCircularKernel,Line_arc_6,Circular_arc_6>  BBCircVariantTraits;

  bench.kernel("CK BBox VarTraits") ;

   bench.Compute<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>(Dxffilename[i]);
  //  bench.Compute_no_dxf<BBCircularKernel,BBCircVariantTraits,BBCircVarContainer>();

 #endif

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

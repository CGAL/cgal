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
{
    char* Dxffilename[]={"myFirst.dxf","minimask0.dxf","minimask1.dxf","mask0.dxf","smallpainttrack.dxf","mask0_25.dxf","mask0_5.dxf","mask1.dxf","cad_l2.dxf","cad_l1.dxf"/*,"CIOnZDraw.dxf","che_mod1.dxf","elekonta.dxf","painttrack.dxf","netlist_signal_1.dxf","51.dxf"*/};
    char* Htmlfilename;
    char* Texfilename;
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
                //Dxffilename[i]="CIOnZDraw.dxf";

                //Dxffilename[]={"myFirst.dxf","minimask0.dxf","minimask1.dxf","mask0.dxf"};
//                 Dxffilename[]={/*"myFirst.dxf","minimask0.dxf","minimask1.dxf","mask0.dxf","smallpainttrack.dxf","mask0_25.dxf","mask0_5.dxf","cad_l2.dxf","cad_l1.dxf",*/"CIOnZDraw.dxf","che_mod1.dxf","elekonta.dxf","painttrack.dxf","netlist_signal_1.dxf","51.dxf"};
                Htmlfilename="benchmarks.html";
                Texfilename="benchmarks.tex";
                std::cout<< "Default *.dxf file is  :" << Dxffilename[i]<< std::endl;
                std::cout<< "Default *.html file is :" << Htmlfilename<< std::endl;
                std::cout<< "Default *.tex file is  :" << Texfilename<< std::endl;
 }




Bench bench(Htmlfilename,Texfilename,Dxffilename[i],true);

 for(i=0;i<15;i++){

 std::ifstream fin;
 fin.open (Dxffilename[i]);
 if (!fin.is_open())
  {
    std::cout<<"file "<< Dxffilename[i] << " is not found"<<std::endl;
    fin.close();
    return 0;

  }
 fin.close();
/*-------------------------------------------------------------------------------------------------------------------------
                                                  !!!!!!!!!!!Circular_Kernel!!!!!!!!!!!!!!!!!!

  -------------------------------------------------------------------------------------------------------------------------*/

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1>     CircularKernel;
/*
  #ifndef CGAL_CURVED_KERNEL_DEBUG
   typedef CGAL::Circular_arc_traits<CircularKernel>                   CircularK_CA_Traits;
  #else
   typedef  CGAL::Circular_arc_traits<CircularKernel>                  Traits0;
   typedef  CGAL::Circular_arc_traits_tracer<Traits0>                  CircularK_CA_Traits;
  #endif

  typedef CircularKernel::Circular_arc_2                                          CircularKArc;
  typedef std::vector<CircularKArc>                                      CircularKArcContainer;
  bench.kernel("Circular kernel  Circular arc traits");

  bench.Compute_no_dxf<CircularKernel,CircularK_CA_Traits,CircularKArcContainer>();*/


  typedef CircularKernel::Circular_arc_2                                  Circular_arc_2;
  typedef CircularKernel::Line_arc_2                                      Line_arc_2;
  typedef CGAL::Variant_traits<CircularKernel,Line_arc_2,Circular_arc_2>  CircularK_Variant_Traits;

  typedef std::variant< Circular_arc_2, Line_arc_2 >        CircularKVarArc;
  typedef std::vector<CircularKVarArc>                        CircularKVarArcContainer;

  bench.kernel("Circular kernel  Variant traits");

//   bench.Compute<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(Dxffilename[i]);
  bench.Compute_dxf<CircularKernel,CircularK_Variant_Traits,CircularKVarArcContainer>(Dxffilename[i]);


  if (i+1<15)
  {
 if (strcmp(Dxffilename[i+1],""))
         {
         bench.newDxfFilename(Dxffilename[i+1]);
         }
 else
         {
        std::cout << "that's all" << Dxffilename[i+1] << std::endl;
        break;
         }
 }
 }

bench.infotable();

  return 0;
}

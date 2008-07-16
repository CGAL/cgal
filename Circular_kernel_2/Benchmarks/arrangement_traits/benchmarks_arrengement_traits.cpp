#include <CGAL/basic.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Algebraic_kernel_2_2.h>

#include <CGAL/intersections.h>

#include <CGAL/Circular_kernel.h>
#include <CGAL/Circular_arc_traits.h>
#include <CGAL/Circular_arc_traits_tracer.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>


#include <CGAL/Random.h>
#include "benchmark.h"

#include <fstream>
#include "IntStr.h"
int main(int argc, char* argv[])
{   char* Dxffilename; 
    char* Htmlfilename;
    char* Texfilename;
    char exten[4];
 

 
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
		Dxffilename = argv[1]; 
		std::cout<< "File "<< Dxffilename << " is correct."<<std::endl;
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
 		Dxffilename="cad_l2.dxf";
		Htmlfilename="benchmarks.html";
		Texfilename="benchmarks.tex";
		std::cout<< "Default *.dxf file is  :" << Dxffilename<< std::endl;
		std::cout<< "Default *.html file is :" << Htmlfilename<< std::endl;
		std::cout<< "Default *.tex file is  :" << Texfilename<< std::endl;
 }
 
std::ifstream fin;
fin.open (Dxffilename);
if (!fin.is_open())
  {
    cout<<"file "<< Dxffilename << " is not found"<<std::endl;
    fin.close();
    return 0;
    
  }
fin.close();
Bench bench(Htmlfilename,Texfilename,Dxffilename);
/*-------------------------------------------------------------------------------------------------------------------------
  						!!!!!!!!!!!Arrangement_traits!!!!!!!!!!!!!!!!!!
						
  -------------------------------------------------------------------------------------------------------------------------*/

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
                                 Alg_kernel,
                                 Nt_traits>           Traits_2;

typedef Traits_2::Curve_2                             Conic_arc_2;
typedef Traits_2::Polygon_2                             Polygon;
bench.kernel("Arrangement_traits");

typedef  std::vector<Conic_arc_2>                     ArcContainer;	
  
bench.Compute_no_dxf<Rat_kernel,Traits_2,ArcContainer>(); 


bench.infotable(); 

  return 0;
}


  


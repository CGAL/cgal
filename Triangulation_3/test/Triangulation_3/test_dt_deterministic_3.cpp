#define CGAL_HILBERT_SORT_WITH_MEDIAN_POLICY_CROSS_PLATFORM_BEHAVIOR 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>


#define GRID_SIZE 10

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;


int main()
{
  std::vector <Kernel::Point_3> points;
  points.reserve(GRID_SIZE*GRID_SIZE*GRID_SIZE);
  
  for (int i=0;i<GRID_SIZE;++i)
    for (int j=0;j<GRID_SIZE;++j)
      for (int k=0;k<GRID_SIZE;++k)
        points.push_back(Kernel::Point_3(i,j,k));
  
  DT3 dt3(points.begin(),points.end());
  
  
  std::stringstream buffer;
  buffer << dt3;
  //~ std::ofstream out ("test_dt_deterministic_3.in");
  //~ out << dt3;
  
  
  //reading the result from a file
  std::ifstream file ("test_dt_deterministic_3.in");
  std::string original,computed;
  if (file)
  {
    while (file && file >> original ){
      buffer >> computed;
      
      if ( original!=computed ){
        std::cout <<"Error: triangulations are differents"<< std::endl;
        std::cout << "|" << original <<"| vs |"<< computed << "|"<< std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  
  std::cout <<"Triangulations are identical"<< std::endl;
  return 0;
}

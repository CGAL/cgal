#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstring>

typedef CGAL::Linear_cell_complex<2,2> LCC_2;
typedef LCC_2::Dart_handle             Dart_handle;
typedef LCC_2::Point                   Point;

int main(int narg, char** argv)
{
  if ( narg!=2 || !strcmp(argv[1],"-h") || !strcmp(argv[1],"-?") )
  {
    std::cout<<"Usage: a.out filename"<<std::endl
             <<"  with filename the name of the file containing a 2D "
             <<"plane graph."<<std::endl<<std::endl
             <<"File must be in text mode, respecting the following format:"
             <<std::endl
             <<"***********************************************"<<std::endl
             <<"nbvertices nbedges"<<std::endl
             <<"x y //coordinates, one pair for each vertex"<<std::endl
             <<"..."<<std::endl
             <<"i j //edge betwen vertices number i and j,"
             <<" one pair for each edge"<<std::endl
             <<"..."<<std::endl
             <<"***********************************************"<<std::endl
             <<std::endl;
    return EXIT_FAILURE;
  }
  
  LCC_2 lcc;

  std::ifstream is(argv[1]);
  std::cout<<"Import plane graph from "<<argv[1]<<std::endl;
  CGAL::import_from_plane_graph(lcc, is);
  
  // Display the lcc characteristics.
  std::cout<<"LCC characteristics:"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) 
    << ", valid=" << lcc.is_valid() << std::endl;
  
  return EXIT_SUCCESS;
}


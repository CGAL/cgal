#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>

/* If you want to use a viewer, you can use one of the following file
 * depending if you use vtk or qglviewer. */
//#include "linear_cell_complex_3_viewer_qt.h"
//#include "linear_cell_complex_3_viewer_vtk.h"

typedef CGAL::Linear_cell_complex<3> LCC_3;
typedef LCC_3::Dart_handle           Dart_handle;
typedef LCC_3::Point                 Point;

typedef CGAL::Delaunay_triangulation_3<LCC_3::Traits> Triangulation;

int main(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cout<<"Usage : voronoi_3 filename"<<std::endl   
             <<"   filename being a fine containing 3D points used to "
             <<" compute the Delaunay_triangulation_3."<<std::endl;
    return EXIT_FAILURE;
  }

  // 1) Compute the Delaunay_triangulation_3.
  Triangulation T;

  std::ifstream iFile(argv[1], std::ios::in);
  if (!iFile)
  {
    std::cout << "Problem reading file " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }
  
  std::istream_iterator<Point> begin(iFile), end;
  T.insert(begin, end);
  assert(T.is_valid(false));
 
  // 2) Convert the triangulation into a 3D lcc.
  LCC_3 lcc;
  CGAL::import_from_triangulation_3<LCC_3, Triangulation>(lcc, T);

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  /* The following line is possible only if you use
   *  cgal_lcc_viewer_qt or cgal_lcc_viewer_vtk */
  //  display_lcc(lcc); 

  // 3) Compute the dual lcc.
  LCC_3 dual_lcc;
  CGAL::dual<LCC_3>(lcc,dual_lcc);
  // Here, dual_lcc is the 3D Voronoi diagram.
  
  // 4) Display the dual_lcc characteristics.
  std::cout<<"Voronoi subdvision :"<<std::endl<<"  ";
  dual_lcc.display_characteristics(std::cout) << ", valid=" 
                                              << dual_lcc.is_valid()
                                              << std::endl;

  return EXIT_SUCCESS;
}


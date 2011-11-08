#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>

/* If you want to use a viewer, you can use one of the following file
 * depending if you use vtk or qglviewer. */
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#else 
#ifdef CGAL_LCC_USE_VTK
#include "linear_cell_complex_3_viewer_vtk.h"
#endif
#endif

typedef CGAL::Linear_cell_complex<3> LCC_3;
typedef LCC_3::Dart_handle           Dart_handle;
typedef LCC_3::Point                 Point;

typedef CGAL::Delaunay_triangulation_3<LCC_3::Traits> Triangulation;

// Function used to display the voronoi diagram.
void display_voronoi(LCC_3& alcc, Dart_handle adart)
{
  // We remove all the volumes containing one dart of the infinite volume
  std::stack<Dart_handle> toremove;
  int mark_toremove=alcc.get_new_mark();

  // We cannot view the infinite volume since it does not have
  // a correct geometry. For that we have to remove the infinite volume.
  toremove.push(adart);
  CGAL::mark_cell<LCC_3,3>(alcc, adart, mark_toremove);
 
  // Plus all the volumes sharing a face with it.
  for (LCC_3::Dart_of_cell_range<3>::iterator
         it=alcc.darts_of_cell<3>(adart).begin(),
         itend=alcc.darts_of_cell<3>(adart).end(); it!=itend; ++it)
  {
    if ( !alcc.is_marked(it->beta(3), mark_toremove) )
    {
      CGAL::mark_cell<LCC_3,3>(alcc, it->beta(3), mark_toremove);
      toremove.push(it->beta(3));
    }    
  }
  
  while( !toremove.empty() )
  {
    CGAL::remove_cell<LCC_3, 3>(alcc, toremove.top());
    toremove.pop();
  }

#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(alcc);
#endif // CGAL_LCC_USE_VIEWER
}

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
  Dart_handle dh=
    CGAL::import_from_triangulation_3<LCC_3, Triangulation>(lcc, T);

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  // 3) Compute the dual lcc.
  LCC_3 dual_lcc;
  Dart_handle ddh=CGAL::dual<LCC_3>(lcc, dual_lcc, dh);
  // Here, dual_lcc is the 3D Voronoi diagram.
  
  // 4) Display the dual_lcc characteristics.
  std::cout<<"Voronoi subdvision :"<<std::endl<<"  ";
  dual_lcc.display_characteristics(std::cout) << ", valid=" 
                                              << dual_lcc.is_valid()
                                              << std::endl;

#ifdef CGAL_LCC_USE_VIEWER
  display_voronoi(dual_lcc, ddh);
#endif // CGAL_LCC_USE_VIEWER

  return EXIT_SUCCESS;
}


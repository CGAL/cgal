#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iostream>
#include <fstream>

/* If you want to use a viewer, you can use one of the following file
 * depending if you use vtk or qglviewer. */
#ifdef CGAL_LCC_USE_QT
//#include "linear_cell_complex_3_viewer_qt.h"
#endif

typedef CGAL::Linear_cell_complex<2> LCC_2;
typedef LCC_2::Dart_handle           Dart_handle;
typedef LCC_2::Point                 Point;

typedef CGAL::Delaunay_triangulation_2<LCC_2::Traits> Triangulation;

// Function used to display the voronoi diagram.
void display_voronoi(LCC_2& alcc, Dart_handle adart)
{
  // We remove all the faces containing one dart of the infinite faces
  std::stack<Dart_handle> toremove;
  int mark_toremove=alcc.get_new_mark();

  // We cannot view the infinite face since it does not have
  // a correct geometry. For that we have to remove the infinite face.
  toremove.push(adart);
  CGAL::mark_cell<LCC_2,2>(alcc, adart, mark_toremove);
 
  // Plus all the faces sharing an edge with it.
  for (LCC_2::Dart_of_cell_range<2>::iterator
         it=alcc.darts_of_cell<2>(adart).begin(),
         itend=alcc.darts_of_cell<2>(adart).end(); it!=itend; ++it)
  {
    if ( !alcc.is_marked(it->beta(2), mark_toremove) )
    {
      CGAL::mark_cell<LCC_2,2>(alcc, it->beta(2), mark_toremove);
      toremove.push(it->beta(2));
    }    
  }
  
  while( !toremove.empty() )
  {
    CGAL::remove_cell<LCC_2, 2>(alcc, toremove.top());
    toremove.pop();
  }

#ifdef CGAL_LCC_USE_VIEWER
  //  display_lcc(alcc);
#endif // CGAL_LCC_USE_VIEWER
}

int main(int narg, char** argv)
{
  if (narg>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?")) )
  {
    std::cout<<"Usage : voronoi_2 filename"<<std::endl   
             <<"   filename being a fine containing 3D points used to "
             <<" compute the Delaunay_triangulation_3."<<std::endl;
    return EXIT_FAILURE;
  }

  std::string filename;
  if ( narg==1 )
  {
    filename=std::string("data/points");
    std::cout<<"No filename given: use data/points by default."<<std::endl;
  }
  else
    filename=std::string(argv[1]);
  
  // 1) Compute the Delaunay_triangulation_2.
  Triangulation T;

  std::ifstream iFile(filename.c_str());
  if (!iFile)
  {
    std::cout << "Problem reading file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  
  std::istream_iterator<Point> begin(iFile), end;
  T.insert(begin, end);
  assert(T.is_valid(false));
 
  // 2) Convert the triangulation into a 2D lcc.
  LCC_2 lcc;
  Dart_handle dh=
    CGAL::import_from_triangulation_2<LCC_2, Triangulation>(lcc, T);

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  // 3) Compute the dual lcc.
  LCC_2 dual_lcc;
  Dart_handle ddh=CGAL::dual<LCC_2>(lcc, dual_lcc, dh);
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


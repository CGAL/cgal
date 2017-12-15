#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_functionalities.h>
#include <vector>

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include "linear_cell_complex_3_viewer_qt.h"
#endif

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<2,3> LCC_3_gmap;

int main(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cout<<"Usage: "<<argv[0]<<" filename.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, argv[1]))
  {
    std::cout<<"PROBLEM reading gile "<<argv[1]<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  CGAL::Combinatorial_map_tools<LCC_3_cmap> cmt(lcc);
  
  cmt.surface_simplification_in_one_vertex();
  std::cout<<"All non loop contracted: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  cmt.initialize_faces();
  cmt.surface_simplification_in_one_face();
  std::cout<<"All faces merges: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;
  
  cmt.surface_quadrangulate();
  std::cout<<"After quadrangulation: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;
  
#ifdef CGAL_USE_BASIC_VIEWER
  //  display_lcc(lcc);
#endif // CGAL_USE_BASIC_VIEWER

  return EXIT_SUCCESS;
}

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Surface_mesh_curve_topology.h>

#include <CGAL/Path_generators.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_lcc_with_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;

///////////////////////////////////////////////////////////////////////////////
int main()
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "data/double-torus-example.off"))
  {
    std::cout<<"ERROR reading file data/double-torus-example.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  CGAL::Surface_mesh_curve_topology<LCC_3_cmap> smct(lcc);

  CGAL::Path_on_surface<LCC_3_cmap> p1(lcc); // A first path
  p1.push_back_by_index(56); // Its starting dart
  for (int i=0; i<3; ++i)
  { p1.extend_positive_turn(2); } // Extend the path
  
  CGAL::Path_on_surface<LCC_3_cmap> p2(lcc); // A second path
  p2.push_back_by_index(202);  // Its starting dart
  for (int i=0; i<3; ++i)
  { p2.extend_negative_turn(2); } // Extend the path
  
  CGAL::Path_on_surface<LCC_3_cmap> p3(lcc); // A third path
  p3.push_back_by_index(411); // Its starting dart
  p3.extend_positive_turn(1); // Extend the path
  for (int i=0; i<3; ++i)
  { p3.extend_positive_turn(2); } 
  p3.extend_positive_turn(1);
  
  bool res1=smct.are_base_point_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")<<" base point homotopic with path p2 (green)."<<std::endl;

  bool res2=smct.are_base_point_homotopic(p1, p3);
  std::cout<<"Path p1 (pink) "<<(res2?"IS":"IS NOT")<<" base point homotopic with path p3 (orange)."<<std::endl;

#ifdef CGAL_USE_BASIC_VIEWER
  std::vector<CGAL::Path_on_surface<LCC_3_cmap> > paths;
  paths.push_back(p1);
  paths.push_back(p2);  
  paths.push_back(p3);  
  CGAL::draw(lcc, paths); // Enable only if CGAL was compiled with Qt5
#endif // CGAL_USE_BASIC_VIEWER
  
  return EXIT_SUCCESS;
}

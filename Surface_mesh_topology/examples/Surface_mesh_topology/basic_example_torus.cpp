#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Homotopy_tester.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;

///////////////////////////////////////////////////////////////////////////////
void create_path_1(CGAL::Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index(0); // Its starting dart
  for (int i=0; i<4; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(CGAL::Path_on_surface<LCC_3_cmap>& p)
{ 
  p.push_back_by_index(1); // Its starting dart
  for (int i=0; i<4; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "data/torus_quad.off"))
  {
    std::cout<<"ERROR reading file data/torus_quad.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  CGAL::Homotopy_tester<LCC_3_cmap> smct(lcc);
  CGAL::Path_on_surface<LCC_3_cmap> p1(lcc), p2(lcc);
  create_path_1(p1);
  create_path_2(p2);
  
  bool res1=smct.is_contractible(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")
           <<" contractible."<<std::endl;
  
  bool res2=smct.are_freely_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res2?"IS":"IS NOT")
           <<" homotopic with path p2 (green)."<<std::endl;

#ifdef CGAL_USE_BASIC_VIEWER
  std::vector<CGAL::Path_on_surface<LCC_3_cmap> > paths={p1, p2};
  CGAL::draw(lcc, paths); // Enable only if CGAL was compiled with Qt5
#endif // CGAL_USE_BASIC_VIEWER
  
  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////


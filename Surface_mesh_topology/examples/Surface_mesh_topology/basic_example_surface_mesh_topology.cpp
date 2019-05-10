#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Homotopy_tester.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_lcc_with_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;

///////////////////////////////////////////////////////////////////////////////
void create_path_1(CGAL::Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index(14); // Its starting dart
  for (int i=0; i<7; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}
///////////////////////////////////////////////////////////////////////////////
void create_path_2(CGAL::Path_on_surface<LCC_3_cmap>& p)
{ 
  for (std::size_t index : {202, 206, 335, 317, 322, 69, 62, 414})
  { p.push_back_by_index(index); }
}
///////////////////////////////////////////////////////////////////////////////
void create_path_3(CGAL::Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index(470); // Its starting dart
  for (int i=0; i<13; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}
///////////////////////////////////////////////////////////////////////////////
int main()
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "data/double-torus-example.off"))
  {
    std::cout<<"ERROR reading file data/double-torus-example.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  CGAL::Homotopy_tester<LCC_3_cmap> smct(lcc);
  CGAL::Path_on_surface<LCC_3_cmap> p1(lcc), p2(lcc), p3(lcc);
  create_path_1(p1);
  create_path_2(p2);
  create_path_3(p3);
  
  bool res1=smct.is_contractible(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")
           <<" contractible."<<std::endl;
  
  bool res2=smct.are_freely_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res2?"IS":"IS NOT")
           <<" homotopic with path p2 (green)."<<std::endl;

  bool res3=smct.are_freely_homotopic(p1, p3);
  std::cout<<"Path p1 (pink) "<<(res3?"IS":"IS NOT")
           <<" homotopic with path p3 (orange)."<<std::endl;

#ifdef CGAL_USE_BASIC_VIEWER
  std::vector<CGAL::Path_on_surface<LCC_3_cmap> > paths={p1, p2, p3};
  CGAL::draw(lcc, paths); // Enable only if CGAL was compiled with Qt5
#endif // CGAL_USE_BASIC_VIEWER
  
  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////


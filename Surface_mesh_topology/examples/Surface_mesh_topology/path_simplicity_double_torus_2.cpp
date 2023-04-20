#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
void create_path(Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index(682); // Its starting dart
  for (int i=0; i<11; ++i)
  { p.extend_positive_turn(2); } // Extend the path
  p.extend_positive_turn(3);
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  p.extend_positive_turn(2);
  p.extend_positive_turn(3);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  p.extend_positive_turn(2);
  p.extend_positive_turn(1);
  for (int i=0; i<8; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  for (int i=0; i<4; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  for (int i=0; i<3; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  p.extend_positive_turn(1);
  for (int i=0; i<11; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  bool draw=(argc>1?std::string(argv[1])=="-draw":false);
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, CGAL::data_file_path("meshes/double-torus-example.off").c_str()))
  {
    std::cout<<"ERROR reading file data/double-torus.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  Curves_on_surface_topology<LCC_3_cmap> cst(lcc);
  Path_on_surface<LCC_3_cmap> p(lcc);
  create_path(p);

  bool res=cst.is_homotopic_to_simple_cycle(p);
  std::cout<<"Path p (pink) "<<(res?"IS":"IS NOT")
           <<" simple."<<std::endl;

  if (draw)
  { CGAL::draw(lcc, {p}); }

  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////


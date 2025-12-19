#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
void create_path_1(Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index(438); // Its starting dart
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); } // Extend the path
  p.extend_positive_turn(1);
  for (int i=0; i<7; ++i)
  { p.extend_positive_turn(2); }
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index({14, 15, 391, 392, 395, 227, 223, 313, 318, 326, 82,
                        87, 431, 160, 435, 753, 754, 756, 757, 674, 678, 850,
                        483, 480, 475, 470, 893, 618, 622, 548, 551, 795,
                        797, 806, 637, 634, 638, 872, 521, 376, 180, 424,
                        88, 95, 440, 152, 149, 21});
}

///////////////////////////////////////////////////////////////////////////////
void create_path_3(Path_on_surface<LCC_3_cmap>& p)
{
  p.push_back_by_index(473); // Its starting dart
  p.extend_positive_turn(1); // Extend the path
  p.extend_positive_turn(3);
  for (int i=0; i<7; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  p.extend_positive_turn(1);
  for (int i=0; i<3; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  p.extend_positive_turn(1);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  for (int i=0; i<3; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  p.extend_positive_turn(2);
  p.extend_positive_turn(1);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  for (int i=0; i<3; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  p.extend_positive_turn(1);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  for (int i=0; i<3; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  for (int i=0; i<3; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(1);
  p.extend_positive_turn(3);
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); }
  p.extend_positive_turn(3);
  p.extend_positive_turn(2);
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
  Path_on_surface<LCC_3_cmap> p1(lcc), p2(lcc), p3(lcc);
  create_path_1(p1);
  create_path_2(p2);
  create_path_3(p3);


  bool res1=cst.is_homotopic_to_simple_cycle(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")
           <<" simple."<<std::endl;

  bool res2=cst.is_homotopic_to_simple_cycle(p2);
  std::cout<<"Path p2 (green) "<<(res2?"IS":"IS NOT")
           <<" simple."<<std::endl;

  bool res3=cst.is_homotopic_to_simple_cycle(p3);
  std::cout<<"Path p3 (blue) "<<(res3?"IS":"IS NOT")
           <<" simple."<<std::endl;

  if (draw)
  {
    auto cycles={p1, p2, p3};
    CGAL::draw(lcc, cycles);
  }

  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////


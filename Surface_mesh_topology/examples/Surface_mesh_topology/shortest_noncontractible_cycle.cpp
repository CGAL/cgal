#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>
#include <CGAL/Random.h>
#include <iostream>
#include <cstdlib>

using LCC_3          =CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using Path_on_surface=CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;

double cycle_length(const LCC_3& lcc, const Path_on_surface& cycle)
{ // Compute the length of the given cycle.
  double res=0;
  for (std::size_t i=0; i<cycle.length(); ++i)
  { res+=std::sqrt
        (CGAL::squared_distance(lcc.point(cycle[i]),
                                lcc.point(lcc.other_extremity(cycle[i])))); }
  return res;
}

void display_cycle_info(const LCC_3& lcc, const Path_on_surface& cycle)
{ // Display information about the given cycle.
  if (cycle.is_empty()) { std::cout<<"Empty."<<std::endl; return; }
  std::cout<<"Root: "<<lcc.point(cycle[0])<<"; "
           <<"Number of edges: "<<cycle.length()<<"; "
           <<"Length: "<<cycle_length(lcc, cycle)<<std::endl;
}

int main(int argc, char* argv[])
{
  std::string filename(argc==1?CGAL::data_file_path("meshes/3torus.off"):argv[1]);
  bool draw=(argc<3?false:std::string(argv[2])=="-draw");

  LCC_3 lcc;
  if (!CGAL::load_off(lcc, filename.c_str())) // Load the off file.
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  std::cout<<"File '"<<filename<<"' loaded. Finding shortest non contractible cycle..."<<std::endl;

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3> cst(lcc);
  LCC_3::Dart_const_descriptor root=lcc.dart_descriptor
    (CGAL::get_default_random().get_int(0, static_cast<int>(lcc.number_of_darts()))); // One dart of the mesh

  Path_on_surface cycle1=
      cst.compute_shortest_non_contractible_cycle_with_base_point(root);

  CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<LCC_3> wf(lcc);
  Path_on_surface cycle2=
      cst.compute_shortest_non_contractible_cycle_with_base_point(root, wf);

  std::cout<<"Cycle 1 (pink): "; display_cycle_info(lcc, cycle1);
  std::cout<<"Cycle 2 (green): "; display_cycle_info(lcc, cycle2);
  if (draw) { CGAL::draw(lcc, {cycle1, cycle2}); }

  return EXIT_SUCCESS;
}

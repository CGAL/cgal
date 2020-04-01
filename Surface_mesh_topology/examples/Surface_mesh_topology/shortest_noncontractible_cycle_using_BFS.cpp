#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

using LCC_3          =CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using CST            =CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;
using Path_on_surface=CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;

int main(int argc, char* argv[])
{
  std::cout<<"Program shortest_noncontractible_cycle_using_BFS started."
           <<std::endl;
  std::string filename(argc==1?"data/3torus.off":argv[1]);
  bool draw=(argc<3?false:std::string(argv[2])=="-draw");
  
  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  LCC_3 lcc;
  CGAL::load_off(lcc, inp);
  std::cout<<"File '"<<filename<<"' loaded. Finding the shortest noncontractible cycle..."<<std::endl;
  
  CST cst(lcc);
  LCC_3::Dart_handle root=lcc.darts().begin(); // Fix one vertex of the mesh
  Path_on_surface cycle=
      cst.compute_shortest_noncontractible_cycle_with_basepoint(root);
  if (cycle.length()==0)
  { std::cout<<"  Cannot find such cycle. Stop."<<std::endl; }
  else
  {
    std::cout<<"  Number of edges in cycle: "<<cycle.length()<<std::endl;
    std::cout<<"  Root: "<<lcc.point(root)<<std::endl;
    if (draw) { CGAL::draw(lcc, cycle); }
  }

  return EXIT_SUCCESS;
}

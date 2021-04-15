#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/draw_face_graph_with_paths.h>

using Mesh           =CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3>;
using Path_on_surface=CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;

double cycle_length(const Mesh& mesh, const Path_on_surface& cycle)
{ // Compute the length of the given cycle.
  double res=0;
  for (std::size_t i=0; i<cycle.length(); ++i)
  { res+=std::sqrt
        (CGAL::squared_distance(mesh.point(mesh.vertex(mesh.edge(cycle[i]), 0)),
                                mesh.point(mesh.vertex(mesh.edge(cycle[i]), 1)))); }
  return res;
}

void display_cycle_info(const Mesh& mesh, const Path_on_surface& cycle)
{ // Display information about the given cycle.
  if (cycle.is_empty()) { std::cout<<"Empty."<<std::endl; return; }
  std::cout<<"Root: "<<mesh.point(mesh.vertex(mesh.edge(cycle[0]), 0))<<"; "
           <<"Number of edges: "<<cycle.length()<<"; "
           <<"Length: "<<cycle_length(mesh, cycle)<<std::endl;
}

int main(int argc, char* argv[])
{
  std::string filename(argc==1?"data/3torus.off":argv[1]);
  bool draw=(argc<3?false:(std::string(argv[2])=="-draw"));
  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  Mesh sm;
  inp>>sm;
  std::cout<<"File '"<<filename<<"' loaded. Finding edge-width of the mesh..."<<std::endl;

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh> cst(sm, true);

  Path_on_surface cycle1=cst.compute_edge_width(true);

  CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<Mesh> wf(sm);
  Path_on_surface cycle2=cst.compute_shortest_non_contractible_cycle(wf, true);

  std::cout<<"Cycle 1 (pink): "; display_cycle_info(sm, cycle1);
  std::cout<<"Cycle 2 (green): "; display_cycle_info(sm, cycle2);
  if (draw) { CGAL::draw(sm, {cycle1, cycle2}); }

  return EXIT_SUCCESS;
}

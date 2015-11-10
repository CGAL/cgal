#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/get_border.h>

#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor  halfedge_descriptor;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh)) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  double target_edge_length = 0.04;
  unsigned int nb_iter = 3;

  std::cout << "Split border...";

    std::vector<halfedge_descriptor> border;
    PMP::get_border(mesh, faces(mesh), std::back_inserter(border));
    PMP::split_long_edges(mesh, border, target_edge_length);

  std::cout << "done." << std::endl;

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  PMP::isotropic_remeshing(mesh,
      faces(mesh),
      target_edge_length,
      PMP::parameters::number_of_iterations(nb_iter)
      .protect_constraints(true)//i.e. protect border, here
      );

  std::cout << "Remeshing done." << std::endl;

  return 0;
}

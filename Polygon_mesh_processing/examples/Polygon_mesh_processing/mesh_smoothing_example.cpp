#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Surface_mesh<K::Point_3>                          Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const char* filename = argc > 1 ? argv[1] : "data/mech-holes-shark.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid .off file." << std::endl;
    return EXIT_FAILURE;
  }

  std::set<Mesh::Vertex_index> constrained_vertices;
  for(Mesh::Vertex_index v : vertices(mesh))
  {
    if(is_border(v, mesh))
      constrained_vertices.insert(v);
  }

  std::cout << "Constraining: " << constrained_vertices.size() << " border vertices" << std::endl;

  const unsigned int nb_iterations = 5;
  CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);

  std::cout << "Smoothing... (" << nb_iterations << " iterations)" << std::endl;
  PMP::smooth(mesh, PMP::parameters::number_of_iterations(nb_iterations)
                                    .use_safety_constraints(false) // authorize all moves
                                    .vertex_is_constrained_map(vcmap));

  std::ofstream output("mesh_smoothed.off");
  output << mesh;

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}

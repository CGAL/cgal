#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid .off file." << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int nb_iterations = (argc > 2) ? std::atoi(argv[2]) : 10;
  const double time = (argc > 3) ? std::atof(argv[3]) : 0.0001;

  std::set<Mesh::Vertex_index> constrained_vertices;
  for(Mesh::Vertex_index v : vertices(mesh))
  {
    if(is_border(v, mesh))
      constrained_vertices.insert(v);
  }

  std::cout << "Constraining: " << constrained_vertices.size() << " border vertices" << std::endl;
  CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);

  std::cout << "Smoothing shape... (" << nb_iterations << " iterations)" << std::endl;
  PMP::smooth_shape(mesh, time, PMP::parameters::number_of_iterations(nb_iterations)
                                                .vertex_is_constrained_map(vcmap));

  std::ofstream output("mesh_shape_smoothed.off");
  output.precision(17);
  output << mesh;

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}

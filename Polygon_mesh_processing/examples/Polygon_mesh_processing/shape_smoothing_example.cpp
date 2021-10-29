#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
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

  CGAL::IO::write_polygon_mesh("mesh_shape_smoothed.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}

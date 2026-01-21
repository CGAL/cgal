#define CGAL_PMP_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.04;
  unsigned int nb_iter = (argc > 3) ? std::stoi(std::string(argv[3])) : 1;

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
                           CGAL::parameters::number_of_iterations(nb_iter)
                           .smoothing_algo(PMP::FAIR)
                           .do_project(false)
                           );

  CGAL::IO::write_polygon_mesh("out.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Remeshing done." << std::endl;

  return 0;
}

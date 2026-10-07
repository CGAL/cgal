#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>

#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_file>\n";
    return 1;
  }

  Surface_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(argv[1], mesh) || mesh.is_empty())
  {
    std::cerr << "Cannot read " << argv[1] << "\n";
    return 1;
  }

  if(!CGAL::is_triangle_mesh(mesh))
    PMP::triangulate_faces(mesh);

  PMP::autorefine(mesh, CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));

  const bool self_intersecting =
    PMP::does_self_intersect<CGAL::Parallel_if_available_tag>(mesh);

  std::cout << "{\n"
            << "  \"Nb_output_vertices\": " << num_vertices(mesh) << ",\n"
            << "  \"Nb_output_edges\": " << num_edges(mesh) << ",\n"
            << "  \"Nb_output_faces\": " << num_faces(mesh) << ",\n"
            << "  \"Is_triangle_mesh\": \"" << (CGAL::is_triangle_mesh(mesh) ? "True" : "False") << "\",\n"
            << "  \"Is_valid\": \"" << (CGAL::is_valid_polygon_mesh(mesh) ? "True" : "False") << "\",\n"
            << "  \"Self_intersecting\": \"" << (self_intersecting ? "True" : "False") << "\",\n"
            << "  \"Closed\": \"" << (CGAL::is_closed(mesh) ? "True" : "False") << "\",\n"
            << "  \"Output_bound_a_volume\": \"" << (PMP::does_bound_a_volume(mesh) ? "True" : "False") << "\"\n"
            << "}\n";

  return 0;
}

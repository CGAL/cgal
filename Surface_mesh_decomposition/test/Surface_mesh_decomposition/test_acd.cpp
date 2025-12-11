#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/approximate_convex_decomposition.h>

#include <iostream>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Point = K::Point_3;

using Convex_hull = std::pair<std::vector<Point>, std::vector<std::array<unsigned int, 3> > >;
using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = Mesh::Vertex_index;
using face_descriptor = Mesh::Face_index;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/knot2.off");
  std::cout << filename << std::endl;

  Mesh mesh;

  std::vector<Convex_hull> convex_volumes;

  // Try with empty mesh
  CGAL::approximate_convex_decomposition(mesh, std::back_inserter(convex_volumes),
    CGAL::parameters::maximum_depth(10)
    .volume_error(0.1)
    .maximum_number_of_convex_volumes(9)
    .split_at_concavity(true)
    .maximum_number_of_voxels(1000000)
    .concurrency_tag(CGAL::Parallel_if_available_tag()));

  assert(convex_volumes.size() == 0);

  vertex_descriptor v0 = mesh.add_vertex(Point(0, 0, 0));
  vertex_descriptor v1 = mesh.add_vertex(Point(1, 0, 0));
  vertex_descriptor v2 = mesh.add_vertex(Point(0, 1, 0));
  vertex_descriptor v3 = mesh.add_vertex(Point(1, 0, 0));
  vertex_descriptor v4 = mesh.add_vertex(Point(0, 1, 0));

  face_descriptor fd = mesh.add_face(v0, v1, v2);
  mesh.add_face(v1, v4, v3);

  CGAL::approximate_convex_decomposition(mesh, std::back_inserter(convex_volumes),
    CGAL::parameters::maximum_depth(10)
    .volume_error(0.1)
    .maximum_number_of_convex_volumes(9)
    .split_at_concavity(true)
    .maximum_number_of_voxels(1000000)
    .concurrency_tag(CGAL::Parallel_if_available_tag()));

  for (std::size_t i = 0; i < convex_volumes.size(); i++) {
    const Convex_hull& ch = convex_volumes[i];
    CGAL::IO::write_polygon_soup(std::to_string(i) + ".off", ch.first, ch.second);
  }

  assert(convex_volumes.size() == 1);

  mesh.clear();

  if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  convex_volumes.clear();
  convex_volumes.reserve(9);

  CGAL::approximate_convex_decomposition(mesh, std::back_inserter(convex_volumes),
    CGAL::parameters::maximum_depth(10)
    .volume_error(0.1)
    .maximum_number_of_convex_volumes(9)
    .split_at_concavity(true)
    .maximum_number_of_voxels(1000000)
    .concurrency_tag(CGAL::Parallel_if_available_tag()));

  assert(convex_volumes.size() > 0);
  assert(convex_volumes.size() <= 9);

  return 0;
}

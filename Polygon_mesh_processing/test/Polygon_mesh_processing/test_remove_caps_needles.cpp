#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor     face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  std::cout << "Input mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "Input mesh has self-intersections\n";

  PMP::experimental::remove_almost_degenerate_faces(mesh,
                                                    std::cos(160. / 180 * CGAL_PI),
                                                    4,
                                                    0.14);


  CGAL::IO::write_polygon_mesh("cleaned_mesh.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Output mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "Output mesh has self-intersections\n";

  return 0;
}

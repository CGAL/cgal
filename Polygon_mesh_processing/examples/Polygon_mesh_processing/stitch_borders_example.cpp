
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/full_border_quads.off";
  std::ifstream input(filename);

  Polyhedron mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  std::cout << "Before stitching : " << std::endl;
  std::cout << "\t Number of vertices  :\t" << mesh.size_of_vertices() << std::endl;
  std::cout << "\t Number of halfedges :\t" << mesh.size_of_halfedges() << std::endl;
  std::cout << "\t Number of facets    :\t" << mesh.size_of_facets() << std::endl;

  CGAL::Polygon_mesh_processing::stitch_borders(mesh);

  std::cout << "Stitching done : " << std::endl;
  std::cout << "\t Number of vertices  :\t" << mesh.size_of_vertices() << std::endl;
  std::cout << "\t Number of halfedges :\t" << mesh.size_of_halfedges() << std::endl;
  std::cout << "\t Number of facets    :\t" << mesh.size_of_facets() << std::endl;

  std::ofstream output("mesh_stitched.off");
  output.precision(17);
  output << std::setprecision(17) << mesh;

  return 0;
}

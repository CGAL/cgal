
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

int main()
{
  Surface_mesh mesh;
  std::ifstream input("data/full_border_quads.off");
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  std::cout << "Before stitching : " << std::endl;
  std::cout << "\t Number of vertices :\t" << mesh.number_of_vertices() << std::endl;
  std::cout << "\t Number of edges :\t" << mesh.number_of_edges() << std::endl;
  std::cout << "\t Number of facets :\t" << mesh.number_of_faces() << std::endl;

  CGAL::Polygon_mesh_processing::stitch_borders(mesh);

  std::cout << "Stitching done : " << std::endl;
  std::cout << "\t Number of vertices :\t" << mesh.number_of_vertices() << std::endl;
  std::cout << "\t Number of edges :\t" << mesh.number_of_edges() << std::endl;
  std::cout << "\t Number of facets :\t" << mesh.number_of_faces() << std::endl;

  return 0;
}

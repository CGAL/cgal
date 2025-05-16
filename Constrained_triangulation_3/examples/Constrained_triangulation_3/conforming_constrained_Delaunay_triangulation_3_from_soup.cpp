#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <CGAL/draw_constrained_triangulation_3.h>

#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

int main(int argc, char* argv[])
{
  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t>> polygons;

  auto filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cubes.off");
  std::ifstream in(filename);
  if(!in || !CGAL::IO::read_OFF(in, points, polygons)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Read " << points.size() << " vertices and "
            << polygons.size() << " polygons" << std::endl;

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);

  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}

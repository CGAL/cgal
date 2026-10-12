#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/draw_constrained_triangulation_3.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <vector>
#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

int main(int argc, char* argv[])
{
  auto filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cubes.off");

  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t>> polygons;
  if(!CGAL::IO::read_polygon_soup(filename, points, polygons)) {
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

  // Collect constrained facets per polygon
  std::vector<std::size_t> constrained_facets(polygons.size());
  for(auto facet : ccdt.constrained_facets())
  {
    int i = ccdt.face_constraint_index(facet);
    ++constrained_facets[i];
  }
  auto it = std::max_element(constrained_facets.begin(), constrained_facets.end());

  std::cout << "The polygon with the most constrained facets has index "
            << (it - constrained_facets.begin()) << " and " << *it << " facets.\n";

  std::ofstream ofs(argc > 2 ? argv[2] : "out.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, ccdt);

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}

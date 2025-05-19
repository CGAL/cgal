#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "CGAL/Polygon_mesh_processing/polygon_soup_self_intersections.h"
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/draw_constrained_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const auto filename = (argc > 1) ? argv[1]
                       : CGAL::data_file_path("meshes/spheres_intersecting.off");

  CGAL::Surface_mesh<K::Point_3> mesh;
  std::ifstream in(filename);
  if(!in || !(in >> mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Number of facets in " << filename << ": "
    << mesh.number_of_faces() << "\n";

  if(PMP::does_self_intersect(mesh))
  {
    std::cout << "Mesh self-intersects, performing autorefine...\n";

    // use a polygon soup as container as the output will most likely be non-manifold
    std::vector<K::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);
    PMP::autorefine_triangle_soup(points, polygons);
    std::cout << "Number of input triangles after autorefine: "
              << polygons.size() << "\n";

    if(PMP::does_polygon_soup_self_intersect(points, polygons)) {
      std::cerr << "Error: mesh still self-intersects after autorefine\n";
      return EXIT_FAILURE;
    }

    auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(points, polygons);

    std::cout << "Number of constrained facets in the CDT: "
              << ccdt.number_of_constrained_facets() << '\n';

    CGAL::draw(ccdt);
  }
  else
  {
    auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh);

    std::cout << "Number of constrained facets in the CDT: "
              << ccdt.number_of_constrained_facets() << '\n';

    CGAL::draw(ccdt);
  }

  return EXIT_SUCCESS;
}

// STL includes.
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>

#include <CGAL/Shape_detection/Region_growing/internal/Polyline_graph.h>
#include "../../examples/Shape_detection/include/utils.h"

namespace SD  = CGAL::Shape_detection;
using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = typename Kernel::Point_3;

using Surface_mesh   = CGAL::Surface_mesh<Point_3>;
using Polyline_graph = CGAL::Shape_detection::internal::Polyline_graph_points<Surface_mesh>;

using Vertex_to_point_map = typename Polyline_graph::Vertex_to_point_map;
using Point_range         = typename Polyline_graph::Point_range;
using Point_map           = typename Polyline_graph::Point_map;

using Region_type    = CGAL::Shape_detection::Polyline::Least_squares_line_fit_region<Kernel, Point_range, Point_map>;
using Sorting        = CGAL::Shape_detection::Polyline::Least_squares_line_fit_sorting<Kernel, Point_range, Polyline_graph, Point_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Point_range, Polyline_graph, Region_type, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/am.off");
  CGAL::set_ascii_mode(in);
  assert(in);

  Surface_mesh surface_mesh;
  in >> surface_mesh;
  in.close();
  // std::cout << "- num faces: " << faces(surface_mesh).size() << std::endl;
  // std::cout << "- num vertices: " << vertices(surface_mesh).size() << std::endl;
  assert(faces(surface_mesh).size() == 7320);
  assert(vertices(surface_mesh).size() == 3662);

  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::internal::region_growing_planes(
    surface_mesh, std::back_inserter(regions));
  // std::cout << "- num regions: " << regions.size() << std::endl;
  assert(regions.size() == 9);

  // std::string fullpath = (argc > 2 ? argv[2] : "regions_sm.ply");
  // utils::save_polygon_mesh_regions(surface_mesh, regions, fullpath);

  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, surface_mesh));

  Polyline_graph pgraph(surface_mesh, regions,
  CGAL::parameters::vertex_point_map(vertex_to_point_map));
  const auto& point_range = pgraph.point_range();
  // std::cout << "- num extracted points: " << point_range.size() << std::endl;

  Region_type region_type(
    point_range, CGAL::parameters::point_map(pgraph.point_map()));
  Sorting sorting(point_range, pgraph, CGAL::parameters::point_map(pgraph.point_map()));
  sorting.sort();

  Region_growing region_growing(
    point_range, pgraph, region_type, sorting.seed_map());

  std::vector< std::vector<std::size_t> > subregions;
  region_growing.detect(std::back_inserter(subregions));
  // std::cout << "- num subregions: " << subregions.size() << std::endl;
  assert(subregions.size() == 21);

  // fullpath = (argc > 2 ? argv[2] : "subregions_sm.ply");
  // utils::save_point_regions_3<Kernel, Point_range, Point_map>(
  //   point_range, subregions, fullpath, pgraph.point_map());

  std::cout << "rg_pgraph, epick_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}

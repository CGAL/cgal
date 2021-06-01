#include <CGAL/IO/PLY.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_segment_set.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include "include/utils.h"

// Typedefs.
using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = typename Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Face_range   = typename Surface_mesh::Face_range;
using Edge_range   = typename Surface_mesh::Edge_range;

using One_ring_query = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Surface_mesh>;
using Plane_region   = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Surface_mesh>;
using RG_planes      = CGAL::Shape_detection::Region_growing<Face_range, One_ring_query, Plane_region>;

using Face_to_region_map = typename Plane_region::Face_to_region_map;
using Polyline_graph     = CGAL::Shape_detection::Polygon_mesh::Polyline_graph<Kernel, Surface_mesh, Face_to_region_map>;
using Segment_range      = typename Polyline_graph::Segment_range;
using Segment_map        = typename Polyline_graph::Segment_map;

using Line_region  = CGAL::Shape_detection::Segment_set::Least_squares_line_fit_region<Kernel, Segment_range, Segment_map>;
using Line_sorting = CGAL::Shape_detection::Segment_set::Least_squares_line_fit_sorting<Kernel, Segment_range, Polyline_graph, Segment_map>;
using RG_lines     = CGAL::Shape_detection::Region_growing<Segment_range, Polyline_graph, Line_region, typename Line_sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  const std::string filename = is_default_input ? "data/am.off" : argv[1];

  Surface_mesh surface_mesh;
  if (!CGAL::IO::read_polygon_mesh(filename, surface_mesh, CGAL::parameters::all_default())) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }
  const Face_range face_range = faces(surface_mesh);
  const Edge_range edge_range = edges(surface_mesh);
  std::cout << "* number of input faces: " << face_range.size() << std::endl;
  std::cout << "* number of input edges: " << edge_range.size() << std::endl;
  assert(is_default_input && face_range.size() == 7320);
  assert(is_default_input && edge_range.size() == 10980);

  // Find planar regions.
  One_ring_query one_ring_query(surface_mesh);
  Plane_region plane_region(surface_mesh, CGAL::parameters::all_default());
  RG_planes rg_planes(face_range, one_ring_query, plane_region);

  std::vector< std::vector<std::size_t> > regions;
  rg_planes.detect(std::back_inserter(regions));
  std::cout << "* number of found planar regions: " << regions.size() << std::endl;
  assert(is_default_input && regions.size() == 9);

  std::string fullpath = (argc > 2 ? argv[2] : "regions_sm.ply");
  utils::save_polygon_mesh_regions(surface_mesh, regions, fullpath);

  // Find linear regions.
  const Face_to_region_map face_to_region_map(face_range, regions);
  Polyline_graph pgraph(surface_mesh, CGAL::parameters::face_index_map(face_to_region_map));
  const auto& segment_range = pgraph.segment_range();
  std::cout << "* number of extracted segments: " << segment_range.size() << std::endl;

  Line_region line_region(
    segment_range, CGAL::parameters::segment_map(pgraph.segment_map()));
  Line_sorting line_sorting(
    segment_range, pgraph, CGAL::parameters::segment_map(pgraph.segment_map()));
  line_sorting.sort();

  RG_lines rg_lines(
    segment_range, pgraph, line_region, line_sorting.seed_map());

  std::vector< std::vector<std::size_t> > subregions;
  rg_lines.detect(std::back_inserter(subregions));
  std::cout << "* number of found linear regions: " << subregions.size() << std::endl;
  assert(is_default_input && subregions.size() == 21);

  fullpath = (argc > 2 ? argv[2] : "subregions_sm.ply");
  utils::save_segment_regions_3<Kernel, Segment_range, Segment_map>(
    segment_range, subregions, fullpath, pgraph.segment_map());

  return EXIT_SUCCESS;
}

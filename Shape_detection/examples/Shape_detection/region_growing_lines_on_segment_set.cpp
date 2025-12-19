#include <CGAL/IO/PLY.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Segment_set.h>
#include <CGAL/Shape_detection/Region_growing/Polygon_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include "include/utils.h"

// Typedefs.
using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = typename Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Face_range   = typename Surface_mesh::Face_range;
using Edge_range   = typename Surface_mesh::Edge_range;

using One_ring_query = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Surface_mesh>;
using Plane_region   = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Surface_mesh>;
using RG_planes      = CGAL::Shape_detection::Region_growing<One_ring_query, Plane_region>;

using Polyline_graph     = CGAL::Shape_detection::Polygon_mesh::Polyline_graph<Surface_mesh>;
using Segment_range      = typename Polyline_graph::Segment_range;
using Segment_map        = typename Polyline_graph::Segment_map;

using Line_region  = CGAL::Shape_detection::Segment_set::Least_squares_line_fit_region<Kernel, Surface_mesh::Edge_index, Segment_map>;
using Line_sorting = CGAL::Shape_detection::Segment_set::Least_squares_line_fit_sorting<Kernel, Surface_mesh::Edge_index, Polyline_graph, Segment_map>;
using RG_lines     = CGAL::Shape_detection::Region_growing<Polyline_graph, Line_region>;

int main(int argc, char *argv[]) {

  // Load data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  const std::string filename = is_default_input ? CGAL::data_file_path("meshes/am.off") : argv[1];

  Surface_mesh surface_mesh;
  if (!CGAL::IO::read_polygon_mesh(filename, surface_mesh)) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }
  const Face_range face_range = faces(surface_mesh);
  const Edge_range edge_range = edges(surface_mesh);
  std::cout << "* number of input faces: " << face_range.size() << std::endl;
  std::cout << "* number of input edges: " << edge_range.size() << std::endl;
  assert(!is_default_input || face_range.size() == 7320);
  assert(!is_default_input || edge_range.size() == 10980);

  // Find planar regions.
  One_ring_query one_ring_query(surface_mesh);
  Plane_region plane_region(surface_mesh);
  RG_planes rg_planes(face_range, one_ring_query, plane_region);

  std::vector<typename RG_planes::Primitive_and_region> regions;
  rg_planes.detect(std::back_inserter(regions));
  std::cout << "* number of found planar regions: " << regions.size() << std::endl;
  assert(!is_default_input || regions.size() == 9);

  std::string fullpath = (argc > 2 ? argv[2] : "regions_sm.ply");
  utils::save_polygon_mesh_regions(surface_mesh, regions, fullpath);

  // Find linear regions.
  Polyline_graph pgraph(surface_mesh, rg_planes.region_map());
  const auto& segment_range = pgraph.segment_range();
  std::cout << "* number of extracted segments: " << segment_range.size() << std::endl;

  Line_region line_region(CGAL::parameters::segment_map(pgraph.segment_map()));
  Line_sorting line_sorting(
    segment_range, pgraph, CGAL::parameters::segment_map(pgraph.segment_map()));
  line_sorting.sort();

  RG_lines rg_lines(
    segment_range, line_sorting.ordered(), pgraph, line_region);

  std::vector<typename RG_lines::Primitive_and_region> subregions;
  rg_lines.detect(std::back_inserter(subregions));
  std::cout << "* number of found linear regions: " << subregions.size() << std::endl;
  assert(!is_default_input || subregions.size() == 21);

  fullpath = (argc > 2 ? argv[2] : "subregions_sm.ply");
  utils::save_segment_regions_3<Kernel, std::vector<typename RG_lines::Primitive_and_region>, Segment_map>(
    subregions, fullpath, pgraph.segment_map());

  return EXIT_SUCCESS;
}

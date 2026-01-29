// #define CGAL_PMP_SNAP_DEBUG_PP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/helper.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <iostream>
#include <fstream>
#include <set>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = EPICK::FT;
using Point_3 = EPICK::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;

void stitch(Mesh& mesh)
{
  std::vector<halfedge_descriptor> boundaries;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(boundaries));
  std::cout << "Snapping " << boundaries.size() << " border(s)" << std::endl;

  for(halfedge_descriptor hd : boundaries)
    PMP::stitch_boundary_cycle(hd, mesh);

  PMP::stitch_borders(mesh);
}

bool test(const std::string& filepath,
          const FT tolerance,
          const std::size_t expected_cc,
          const std::size_t expected_boundaries)
{
  std::cout << "== Test " << filepath << " ==" << std::endl;

  Mesh mesh;
  if (!CGAL::IO::read_polygon_mesh(filepath, mesh)) {
    std::cerr << "Error: cannot open mesh file " << filepath << std::endl;
    assert(false);
    return false;
  }

  auto tol_map = mesh.template add_property_map<vertex_descriptor, FT>("v:tol", 0).first;

  std::vector<halfedge_descriptor> boundaries;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(boundaries));
  std::size_t nb_boundaries = boundaries.size();

  std::size_t nb_cc = PMP::internal::number_of_connected_components(mesh);

  std::cout << "Input: " << filepath << "\n"
            << "  Connected components: " << nb_cc
            << ", Boundaries: " << nb_boundaries << std::endl;

  CGAL::Constant_property_map<vertex_descriptor, FT> def_tol(tolerance);

  // Each CC independently

  for (halfedge_descriptor hd : boundaries) {
    std::vector<halfedge_descriptor> border;
    halfedge_descriptor done = hd;
    do {
      border.push_back(hd);
      hd = next(hd, mesh);
    }
    while(hd != done);

    PMP::internal::assign_tolerance_with_local_edge_length_bound(border, tol_map, def_tol, mesh);

    PMP::internal::snap_non_conformal(border, mesh, tol_map, border, mesh, tol_map, true /*self snapping*/,
                                      CGAL::parameters::default_values(), CGAL::parameters::default_values());

    PMP::stitch_boundary_cycle(hd, mesh);
  }

  // Now all CCs at once

  std::vector<halfedge_descriptor> border_vertices;
  PMP::border_halfedges(mesh, std::back_inserter(border_vertices));

  PMP::internal::assign_tolerance_with_local_edge_length_bound(border_vertices, tol_map, def_tol, mesh);
  PMP::internal::snap_non_conformal(border_vertices, mesh, tol_map, border_vertices, mesh, tol_map, true /*self snapping*/,
                                    CGAL::parameters::default_values(), CGAL::parameters::default_values());
  CGAL_assertion(is_valid_polygon_mesh(mesh));

  PMP::stitch_borders(mesh);

  CGAL::IO::write_polygon_mesh("last.ply", mesh, CGAL::parameters::stream_precision(17));

  // Recount boundaries
  boundaries.clear();
  PMP::extract_boundary_cycles(mesh, std::back_inserter(boundaries));
  nb_boundaries = boundaries.size();

  // Recount connected components
  nb_cc = PMP::internal::number_of_connected_components(mesh);

  std::cout << "Output: " << filepath << "\n"
            << "  Connected components: " << nb_cc << " (expected " << expected_cc
            << "), Boundaries: " << nb_boundaries << " (expected " << expected_boundaries << ")" << std::endl;

  assert(nb_boundaries == expected_boundaries);
  assert(nb_cc == expected_cc);

  return true;
}

int main(int, char**)
{
  // for each line: input, base tolerance, expected CC, expected boundary count

  test("data_snapping/part-0.ply", 0.1, 1, 1);
  test("data_snapping/part-1.ply", 0.1, 1, 2);
  test("data_snapping/part-2.ply", 0.1, 1, 1);
  test("data_snapping/part-3.ply", 0.1, 1, 1);
  test("data_snapping/part-4.ply", 0.1, 1, 1);
  test("data_snapping/part-5.ply", 0.1, 1, 1);
  test("data_snapping/part-6.ply", 0.1, 1, 1);
  test("data_snapping/part-7.ply", 0.1, 1, 1);
  test("data_snapping/part-8.ply", 0.1, 1, 1);
  test("data_snapping/part-9.ply", 0.1, 1, 1);

  std::cout << "OK" << std::endl;

  return EXIT_SUCCESS;
}

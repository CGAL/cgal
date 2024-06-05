// #define CGAL_PMP_STITCHING_DEBUG_PP

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Dynamic_property_map.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <algorithm>
#include <deque>
#include <iostream>
#include <iterator>
#include <fstream>
#include <set>
#include <unordered_map>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Mesh>
void test_stitch_boundary_cycles(const std::string fname,
                                 const std::size_t expected_n)
{
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor     halfedge_descriptor;

  std::cout << "Testing stitch_boundary_cycles(); file: " << fname << "..." << std::flush;

  std::ifstream input(fname);
  Mesh mesh;
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return;
  }

  // stitching a single cycle should work from any border halfedge
  for(halfedge_descriptor h : halfedges(mesh))
  {
    if(!is_border(h, mesh))
      continue;

    std::unordered_map<halfedge_descriptor, halfedge_descriptor> h2h;

    Mesh mesh_cpy;
    CGAL::copy_face_graph(mesh, mesh_cpy,
                          CGAL::parameters::halfedge_to_halfedge_output_iterator(std::inserter(h2h, h2h.end())));

    assert(is_border(h2h.at(h), mesh_cpy));

    PMP::stitch_boundary_cycle(h2h.at(h), mesh_cpy);
    assert(is_valid_polygon_mesh(mesh_cpy));
  }

  std::size_t res = PMP::stitch_boundary_cycles(mesh);
  std::cout << "res: " << res << " (expected: " << expected_n << ")" << std::endl;

  assert(res == expected_n);
  assert(is_valid_polygon_mesh(mesh));

  // Just to test the API
  PMP::stitch_boundary_cycles(mesh, params::vertex_point_map(get(CGAL::vertex_point, mesh)));
}

template <typename Mesh>
void test_stitch_boundary_cycles()
{
  test_stitch_boundary_cycles<Mesh>("data_stitching/boundary_cycle.off", 4);
  test_stitch_boundary_cycles<Mesh>("data_stitching/boundary_cycle_2.off", 2);
  test_stitch_boundary_cycles<Mesh>("data_stitching/complex_hole.off", 3);
  test_stitch_boundary_cycles<Mesh>("data_stitching/folded_cycle.off", 2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Mesh>
typename boost::graph_traits<Mesh>::halfedge_descriptor get_border_halfedge(const int edge_id,
                                                                            const Mesh& mesh)
{
  assert(edge_id < static_cast<int>(num_edges(mesh)));

  // id is of the edge because it's easier to
  typename boost::graph_traits<Mesh>::edge_iterator eit = edges(mesh).begin();
  std::advance(eit, edge_id);

  typename boost::graph_traits<Mesh>::halfedge_descriptor h = halfedge(*eit, mesh);
  if(!is_border(h, mesh))
    h = opposite(h, mesh);

  assert(is_border(h, mesh));
  return h;
}

template <typename Mesh>
void test_stitch_borders(const std::string fname,
                         const std::size_t expected_n,
                         const bool per_cc = false,
                         std::set<int> unconstrained_edge_ids = { }, // constrained edges must appear in the output
                         std::set<int> cycle_rep_ids = { }) // restrict stitching to cycles containing these edges
{
  std::cout << "Testing stitch_borders(); file: " << fname << "..." << std::flush;

  std::ifstream input(fname);
  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Error: can not read file.";
    return;
  }

  // Mark some edges as UNCONSTRAINED, the stitched mesh must not contain these edges
  typedef typename boost::graph_traits<Mesh>::edge_descriptor         edge_descriptor;
  typedef CGAL::dynamic_edge_property_t<bool>                         Edge_property_tag;
  typedef typename boost::property_map<Mesh, Edge_property_tag>::type Marked_edges;
  typedef PMP::internal::Halfedges_keeper_with_marked_edge_priority<Marked_edges, Mesh> Keeper;

  Marked_edges marks = get(Edge_property_tag(), mesh);
  int eid = 0;
  for(edge_descriptor e : edges(mesh))
    put(marks, e, (unconstrained_edge_ids.count(eid++) == 0));

  Keeper kpr(marks, mesh);

  // Restrict to given cycles
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor     halfedge_descriptor;
  std::set<halfedge_descriptor> cycle_reps;
  for(const int id : cycle_rep_ids)
    cycle_reps.insert(get_border_halfedge(id, mesh));

  std::size_t res = -1;
  if(cycle_reps.empty())
  {
    res = PMP::stitch_borders(mesh, params::apply_per_connected_component(per_cc)
                                           .halfedges_keeper(kpr));
  }
  else
  {
    res = PMP::stitch_borders(cycle_reps, mesh, params::apply_per_connected_component(per_cc)
                                                       .halfedges_keeper(kpr));
  }

  std::cout << "res: " << res << " (expected: " << expected_n << ")" << std::endl;

  for(edge_descriptor e : edges(mesh)) {
    assert(get(marks, e)); // must not be marked as unconstrained
  }

  assert(res == expected_n);
  assert(is_valid_polygon_mesh(mesh));

  // Just to test the API
  Mesh dummy_mesh;
  std::deque<halfedge_descriptor> empty_deque;
  PMP::stitch_borders(empty_deque, dummy_mesh);
  PMP::stitch_borders(dummy_mesh, params::apply_per_connected_component(true));
  PMP::stitch_borders(dummy_mesh);
  std::deque<std::pair<halfedge_descriptor,halfedge_descriptor>> empty_deque_pair;
  PMP::stitch_borders(dummy_mesh, empty_deque_pair);
}

template <typename Mesh>
void test_stitch_borders()
{
  test_stitch_borders<Mesh>("data_stitching/complex_hole.off", 3, false, {}, {83});
  test_stitch_borders<Mesh>("data_stitching/pinched.off", 2, false, {130, 94});
  test_stitch_borders<Mesh>("data_stitching/pinched.off", 2, false, {130, 94}, {94});
  test_stitch_borders<Mesh>("data_stitching/pinched.off", 0, false, {}, {140}); // outer border, nothing to stitch
  test_stitch_borders<Mesh>("data_stitching/full_border.off", 4);
  test_stitch_borders<Mesh>(CGAL::data_file_path("meshes/quads_to_stitch.off"), 4);
  test_stitch_borders<Mesh>("data_stitching/half_border.off", 2, false, {23, 15});
  test_stitch_borders<Mesh>("data_stitching/incidence_3.off", 3);
  test_stitch_borders<Mesh>("data_stitching/incoherent_patch_orientation.off", 1);
  test_stitch_borders<Mesh>("data_stitching/mid_border.off", 2);
  test_stitch_borders<Mesh>("data_stitching/multiple_incidence.off", 10, false, {3, 5, 16, 24, 31, 45, 53, 65});
  test_stitch_borders<Mesh>("data_stitching/non_stitchable.off", 0);
  test_stitch_borders<Mesh>("data_stitching/non_manifold.off", 0);
  test_stitch_borders<Mesh>("data_stitching/non_manifold2.off", 0);
  test_stitch_borders<Mesh>("data_stitching/two_patches.off", 3);
  test_stitch_borders<Mesh>("data_stitching/nm_cubes.off", 4, true /*per cc*/);
}

template <typename Mesh>
void test_local_stitch_borders()
{
  test_stitch_borders<Mesh>("data_stitching/pinched.off", 2, false, {130, 94});

}

////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Mesh>
void test_degenerate()
{
  std::cout << "Testing degenerate cases" << std::endl;

  typedef typename boost::property_traits<
    typename CGAL::GetVertexPointMap<Mesh>::const_type>::value_type Point;

  Mesh tm;

  CGAL::make_triangle(Point(0,0,0), Point(1,0,0), Point(0,1,0), tm);
  CGAL::make_triangle(Point(0,0,0), Point(1,0,0), Point(0,1,0), tm);
  CGAL::make_triangle(Point(0,0,0), Point(1,0,0), Point(0,1,0), tm);
  CGAL::make_triangle(Point(0,0,0), Point(1,0,0), Point(0,1,0), tm);

  std::size_t res = CGAL::Polygon_mesh_processing::stitch_borders(tm);
  assert(res == 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Mesh>
void test()
{
  test_stitch_boundary_cycles<Mesh>();
  test_stitch_borders<Mesh>();
  test_degenerate<Mesh>();
}

int main()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPICK;

  std::cout << "--- Test EPICK SM" << std::endl;
  test<CGAL::Surface_mesh<EPICK::Point_3> >();

  std::cout << "--- Test EPICK Polyhedron" << std::endl;
  test<CGAL::Polyhedron_3<EPICK> >();

  // Disabled because it takes too long for some test plateforms
//  typedef CGAL::Exact_predicates_exact_constructions_kernel       EPECK;

//  std::cout << "--- Test EPECK SM" << std::endl;
//  test<CGAL::Surface_mesh<EPECK::Point_3> >();

//  std::cout << "--- Test EPECK Polyhedron" << std::endl;
//  test<CGAL::Polyhedron_3<EPECK> >();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}

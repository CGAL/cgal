// Tests that erase_edge() of the topology traits keeps the boundary edges, in
// list and vector mode.  The edges are tagged with their creation index; the
// left and right boundary edges are set to the edges tagged 3 and 2,
// respectively, and the last edge and an interior edge are erased.

#include <cstdlib>

#include "test_utils.h"
#include <CGAL/Arrangement_on_curve_1/Closed_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Identified_topology_traits.h>

template <typename Topology_traits>
void test(const std::string& name, std::size_t& failures) {
  Topology_traits tt;
  auto dm = tt.edge_data_map();
  std::vector<typename Topology_traits::Edge_descriptor> edges{tt.left_boundary_edge()};
  for (int i = 1; i < 5; ++i) edges.push_back(tt.create_edge());
  for (int i = 0; i < 5; ++i) put(dm, edges[i], i);
  tt.set_left_boundary_edge(edges[3]);
  tt.set_right_boundary_edge(edges[2]);
  tt.erase_edge(edges[4]);      // the last edge
  tt.erase_edge(edges[1]);      // an interior edge; shifts the indices in vector mode
  bool ok = (tt.number_of_edges() == 3) && (get(dm, tt.left_boundary_edge()) == 3) &&
    (get(dm, tt.right_boundary_edge()) == 2);
  std::cout << (report(ok, name, "boundary edges", failures) ? "PASS " : "FAIL ") << name << std::endl;
}

int main() {
  std::size_t failures = 0;
  test<Aoc::Identified_topology_traits<Identified_traits, void, int, false>>("list identified", failures);
  test<Aoc::Identified_topology_traits<Identified_traits, void, int, true>>("vector identified", failures);
  test<Aoc::Closed_topology_traits<Closed_traits, void, int, false>>("list closed", failures);
  test<Aoc::Closed_topology_traits<Closed_traits, void, int, true>>("vector closed", failures);
  return (failures == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

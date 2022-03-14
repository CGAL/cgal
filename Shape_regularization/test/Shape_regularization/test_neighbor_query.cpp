#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Segments/Delaunay_neighbor_query_2.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_neighbor_query() {

  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(0, 0), Point_2(4, 0)), // the outer square
    Segment_2(Point_2(4, 0), Point_2(4, 4)),
    Segment_2(Point_2(4, 4), Point_2(0, 4)),
    Segment_2(Point_2(0, 4), Point_2(0, 0)),

    Segment_2(Point_2(1, 1), Point_2(3, 1)), // the inner square
    Segment_2(Point_2(3, 1), Point_2(3, 3)),
    Segment_2(Point_2(3, 3), Point_2(1, 3)),
    Segment_2(Point_2(1, 3), Point_2(1, 1))
  };

  // saver.export_eps_segments(segments, "nq_input", 100);

  std::vector<Indices> groups(2);
  groups[0] = {0, 1, 2, 3}; // external square
  groups[1] = {4, 5, 6, 7}; // internal square

  Indices neighbors;
  NQ neighbor_query(
    segments, CGAL::parameters::default_values());

  // Check unique group.
  Segments edges;
  neighbor_query.get_edges(edges);
  // saver.export_eps_segments(edges, "nq_graph0", 100);
  assert(edges.size() == 17);
  assert(neighbor_query.number_of_groups() == 1);
  assert(neighbor_query.number_of_neighbors() == edges.size() * 2);

  std::size_t total = 0;
  neighbor_query(0, neighbors); assert(neighbors.size() == 5); total += neighbors.size();
  neighbor_query(1, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(2, neighbors); assert(neighbors.size() == 4); total += neighbors.size();
  neighbor_query(3, neighbors); assert(neighbors.size() == 4); total += neighbors.size();
  neighbor_query(4, neighbors); assert(neighbors.size() == 4); total += neighbors.size();
  neighbor_query(5, neighbors); assert(neighbors.size() == 5); total += neighbors.size();
  neighbor_query(6, neighbors); assert(neighbors.size() == 5); total += neighbors.size();
  neighbor_query(7, neighbors); assert(neighbors.size() == 4); total += neighbors.size();
  assert(total == neighbor_query.number_of_neighbors());

  // Check clear.
  neighbor_query.clear();
  neighbor_query.get_edges(edges);
  // saver.export_eps_segments(edges, "nq_graph1", 100);
  assert(edges.size() == 0);
  assert(neighbor_query.number_of_groups() == 0);
  assert(neighbor_query.number_of_neighbors() == edges.size() * 2);

  // Add first group.
  neighbor_query.clear();
  neighbor_query.add_group(groups[0]);
  neighbor_query.get_edges(edges);
  // saver.export_eps_segments(edges, "nq_graph2", 100);
  assert(edges.size() == 5);
  assert(neighbor_query.number_of_groups() == 1);
  assert(neighbor_query.number_of_neighbors() == edges.size() * 2);

  total = 0;
  neighbor_query(0, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(1, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  neighbor_query(2, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(3, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  assert(total == neighbor_query.number_of_neighbors());

  // Add second group.
  neighbor_query.clear();
  neighbor_query.add_group(groups[1]);
  neighbor_query.get_edges(edges);
  // saver.export_eps_segments(edges, "nq_graph3", 100);
  assert(edges.size() == 5);
  assert(neighbor_query.number_of_groups() == 1);
  assert(neighbor_query.number_of_neighbors() == edges.size() * 2);

  total = 0;
  neighbor_query(4, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(5, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  neighbor_query(6, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(7, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  assert(total == neighbor_query.number_of_neighbors());

  // Add groups consequently.
  neighbor_query.clear();
  neighbor_query.add_group(groups[0]);
  neighbor_query.add_group(groups[1]);
  neighbor_query.get_edges(edges);
  // saver.export_eps_segments(edges, "nq_graph4", 100);
  assert(edges.size() == 10);
  assert(neighbor_query.number_of_groups() == 2);
  assert(neighbor_query.number_of_neighbors() == edges.size() * 2);

  total = 0;
  neighbor_query(0, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(1, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  neighbor_query(2, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(3, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  neighbor_query(4, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(5, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  neighbor_query(6, neighbors); assert(neighbors.size() == 3); total += neighbors.size();
  neighbor_query(7, neighbors); assert(neighbors.size() == 2); total += neighbors.size();
  assert(total == neighbor_query.number_of_neighbors());

  // Check list with minimum 2 items.
  neighbor_query.clear();
  const std::list<std::size_t> mini = {0, 1};
  neighbor_query.add_group(mini);
  neighbor_query.get_edges(edges);
  // saver.export_eps_segments(edges, "nq_graph5", 100);
  assert(edges.size() == 1);
  assert(neighbor_query.number_of_groups() == 1);
  assert(neighbor_query.number_of_neighbors() == edges.size() * 2);

  total = 0;
  neighbor_query(0, neighbors); assert(neighbors.size() == 1); total += neighbors.size();
  neighbor_query(1, neighbors); assert(neighbors.size() == 1); total += neighbors.size();
  assert(total == neighbor_query.number_of_neighbors());
}

int main() {
  test_neighbor_query< CGAL::Simple_cartesian<double> >();
  test_neighbor_query< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_neighbor_query< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_neighbor_query: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}

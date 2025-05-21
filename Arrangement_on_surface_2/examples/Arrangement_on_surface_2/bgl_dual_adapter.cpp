//! \file examples/Arrangement_on_surface_2/bgl_dual_adapter.cpp
// Adapting the dual of an arrangement to a BGL graph.

#include <CGAL/config.h>

#include <CGAL/boost/graph/breadth_first_search.h>

#include <boost/graph/visitors.hpp>

#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/graph_traits_dual_arrangement_2.h>
#include <CGAL/Arr_face_index_map.h>

#include "Extended_face_property_map.h"
#include "arr_exact_construction_segments.h"
#include "arr_print.h"

using Dcel = CGAL::Arr_face_extended_dcel<Traits, unsigned int>;
using Ex_arrangement = CGAL::Arrangement_2<Traits, Dcel>;
using Dual_arrangement = CGAL::Dual<Ex_arrangement>;
using Face_index_map = CGAL::Arr_face_index_map<Ex_arrangement>;
using Face_property_map =
  Extended_face_property_map<Ex_arrangement,unsigned int>;

int main() {
  // Construct an arrangement of seven intersecting line segments.
  Point p1(1, 1), p2(1, 4), p3(2, 2), p4(3, 7), p5(4, 4), p6(7, 1), p7(9, 3);
  Ex_arrangement  arr;
  insert(arr, Segment(p1, p6));
  insert(arr, Segment(p1, p4));  insert(arr, Segment(p2, p6));
  insert(arr, Segment(p3, p7));  insert(arr, Segment(p3, p5));
  insert(arr, Segment(p6, p7));  insert(arr, Segment(p4, p7));

  // Create a mapping of the arrangement faces to indices.
  Face_index_map index_map(arr);

  // Perform breadth-first search from the unbounded face, using the event
  // visitor to associate each arrangement face with its discover time.
  int time = -1;
  boost::breadth_first_search(Dual_arrangement(arr), arr.unbounded_face(),
                              boost::vertex_index_map(index_map).visitor
                              (boost::make_bfs_visitor
                               (stamp_times(Face_property_map(), time,
                                            boost::on_discover_vertex()))));

  // Print the discover time of each arrangement face.
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    std::cout << "Discover time " << fit->data() << " for ";
    if (fit != arr.unbounded_face()) {
      std::cout << "face ";
      print_ccb<Ex_arrangement>(fit->outer_ccb());
    }
    else std::cout << "the unbounded face.\n";
  }
  return 0;
}

//! \file examples/Arrangement_on_surface_2/vertical_decomposition.cpp
// Performing vertical decomposition of an arrangement.

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Arr_vertical_decomposition_2.h>

#include "arr_linear.h"

typedef boost::variant<Vertex_const_handle, Halfedge_const_handle,
                         Face_const_handle>              Cell_type;
typedef boost::optional<Cell_type>                       Vert_decomp_type;
typedef std::pair<Vert_decomp_type, Vert_decomp_type>    Vert_decomp_pair;
typedef std::pair<Vertex_const_handle, Vert_decomp_pair> Vert_decomp_entry;
typedef std::list<Vert_decomp_entry>                     Vert_decomp_list;

int main() {
  // Construct the arrangement.
  Arrangement arr;

  insert_non_intersecting_curve(arr, Segment(Point(1, 1), Point(3, 0)));
  insert_non_intersecting_curve(arr, Segment(Point(1, 1), Point(2, 2)));
  insert_non_intersecting_curve(arr, Segment(Point(2, 2), Point(3, 0)));
  insert_non_intersecting_curve(arr, Segment(Point(2, 2), Point(5, 0)));
  insert_non_intersecting_curve(arr, Segment(Point(3, 2), Point(5, 0)));
  insert_non_intersecting_curve(arr, Segment(Point(2, 3), Point(3, 3)));
  insert_non_intersecting_curve(arr, Segment(Point(0, 3), Point(6, 4)));
  insert_non_intersecting_curve(arr, Segment(Point(4, 4), Point(4, 5)));

  // Perform vertical ray-shooting from every vertex and locate the feature
  // that lie below it and the feature that lies above it.
  Vert_decomp_list vd_list;
  CGAL::decompose(arr, std::back_inserter(vd_list));

  // Print the results.
  for (auto vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
    const Vert_decomp_pair& curr = vd_iter->second;
    std::cout << "Vertex (" << vd_iter->first->point() << ") : ";

    std::cout << " feature below: ";
    if (! curr.first) std::cout << "EMPTY";
    else {
      auto* vh = boost::get<Vertex_const_handle>(&*(curr.first));;
      if (vh) std::cout << '(' << (*vh)->point() << ')';
      else {
        auto* hh = boost::get<Halfedge_const_handle>(&*(curr.first));
        CGAL_assertion(hh);
        if (! (*hh)->is_fictitious())
          std::cout << '[' << (*hh)->curve() << ']';
        else std::cout << "NONE";
      }
    }

    std::cout << "   feature above: ";
    if (! curr.second) std::cout << "EMPTY\n";
    else {
      auto* vh = boost::get<Vertex_const_handle>(&*(curr.second));;
      if (vh) std::cout << '(' << (*vh)->point() << ")\n";
      else {
        auto* hh = boost::get<Halfedge_const_handle>(&*(curr.second));
        CGAL_assertion(hh);
        if (! (*hh)->is_fictitious())
          std::cout << '[' << (*hh)->curve() << "]\n";
        else std::cout << "NONE\n";
      }
    }
  }

  return 0;
}

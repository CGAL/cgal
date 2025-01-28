//! \file examples/Arrangement_on_surface_2/consolidated_curve_data.cpp
// Associating a color attribute with segments using the consolidated
// curve-data traits.

#include <CGAL/basic.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>

#include "arr_exact_construction_segments.h"

enum Segment_color {RED, BLUE};

using Data_traits =
  CGAL::Arr_consolidated_curve_data_traits_2<Traits, Segment_color>;
using Colored_segment = Data_traits::Curve_2;
using Colored_arr = CGAL::Arrangement_2<Data_traits>;

int main() {
  Colored_arr arr;

  // Construct an arrangement containing three RED line segments.
  insert(arr, Colored_segment(Segment(Point(-1, -1), Point(1, 3)), RED));
  insert(arr, Colored_segment(Segment(Point(2, 0), Point(3, 3)), RED));
  insert(arr, Colored_segment(Segment(Point(0, 3), Point(2, 5)), RED));

  // Insert three BLUE line segments.
  insert(arr, Colored_segment(Segment(Point(-1, 3), Point(4, 1)), BLUE));
  insert(arr, Colored_segment(Segment(Point(-1, 0), Point(4, 1)), BLUE));
  insert(arr, Colored_segment(Segment(Point(-2, 1), Point(1, 4)), BLUE));

  // Go over all vertices and print just the ones corresponding to intersection
  // points between RED segments and BLUE segments. Skip endpoints of
  // overlapping sections.
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    // Go over the current-vertex incident-halfedges and examine their colors.
    bool       has_red = false, has_blue = false;
    Colored_arr::Halfedge_around_vertex_const_circulator eit, first;
    eit = first = vit->incident_halfedges();
    do {
      // Get the color of the current halfedge.
      if (eit->curve().data().size() == 1) {
        Segment_color color = eit->curve().data().front();
        if (color == RED)       has_red = true;
        else if (color == BLUE) has_blue = true;
      }
    } while (++eit != first);

    // Print the vertex only if incident RED and BLUE edges were found.
    if (has_red && has_blue) {
      std::cout << "Red intersect blue at (" << vit->point() << ")\n";
    }
  }

  // Locate the edges that correspond to a red-blue overlap.
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Go over the incident edges of the current vertex and examine their colors.
    bool has_red{false}, has_blue{false};

    for (auto it = eit->curve().data().begin(); it != eit->curve().data().end();
         ++it)
    {
      if (*it == RED) has_red = true;
      else if (*it == BLUE) has_blue = true;
    }

    // Print the edge only if it corresponds to a red-blue overlap.
    if (has_red && has_blue)
      std::cout << "Red overlap blue at [" << eit->curve() << "]\n";
  }

  return 0;
}

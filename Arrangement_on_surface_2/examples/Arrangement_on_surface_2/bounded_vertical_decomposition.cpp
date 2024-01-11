//! \file examples/Arrangement_on_surface_2/vertical_decomposition.cpp
// Performing vertical decomposition of an arrangement.

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Arr_vertical_decomposition_2.h>

#include "arr_exact_construction_segments.h"

using Object_pair = std::pair<CGAL::Object, CGAL::Object>;
using Vert_decomp_entry = std::pair<Vertex_const_handle, Object_pair>;
using Vert_decomp_list = std::list<Vert_decomp_entry>;

int main() {
  // Construct the arrangement.
  Segment segments[] = {Segment(Point(0, 0), Point(3, 3)),
                        Segment(Point(3, 3), Point(6, 0)),
                        Segment(Point(2, 0), Point(5, 3)),
                        Segment(Point(5, 3), Point(8, 0))};
  Arrangement arr;
  insert(arr, segments, segments + sizeof(segments)/sizeof(Segment));

  // Perform vertical ray-shooting from every vertex and locate the feature
  // that lie below it and the feature that lies above it.
  Vert_decomp_list vd_list;
  CGAL::decompose(arr, std::back_inserter(vd_list));

  // Print the results.
  for (auto vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
    const Object_pair& curr = vd_iter->second;
    std::cout << "Vertex (" << vd_iter->first->point() << ") : ";

    Vertex_const_handle vh;
    Halfedge_const_handle hh;
    Face_const_handle fh;

    std::cout << " feature below: ";
    if (CGAL::assign(hh, curr.first)) std::cout << '[' << hh->curve() << ']';
    else if (CGAL::assign(vh, curr.first))
      std::cout << '(' << vh->point() << ')';
    else if (CGAL::assign(fh, curr.first)) std::cout << "NONE";
    else std::cout << "EMPTY";

    std::cout << "   feature above: ";
    if (CGAL::assign(hh, curr.second))
      std::cout << '[' << hh->curve() << "]\n";
    else if (CGAL::assign(vh, curr.second))
      std::cout << '(' << vh->point() << ")\n";
    else if (CGAL::assign(fh, curr.second)) std::cout << "NONE\n";
    else std::cout << "EMPTY\n";
  }

  return 0;
}

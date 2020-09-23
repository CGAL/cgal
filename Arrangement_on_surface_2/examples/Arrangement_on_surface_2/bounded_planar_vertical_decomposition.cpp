//! \file examples/Arrangement_on_surface_2/vertical_decomposition.cpp
// Performing vertical decomposition of an arrangement.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>               Traits_2;
typedef Traits_2::Point_2                                Point_2;
typedef Traits_2::X_monotone_curve_2                     Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                    Arrangement_2;
typedef Arrangement_2::Vertex_const_handle               Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle             Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle                 Face_const_handle;

typedef boost::variant<Vertex_const_handle, Halfedge_const_handle,
                         Face_const_handle>              Cell_type;
typedef boost::optional<Cell_type>                       Vert_decomp_type;
typedef std::pair<Vert_decomp_type, Vert_decomp_type>    Vert_decomp_pair;
typedef std::pair<Vertex_const_handle, Vert_decomp_pair> Vert_decomp_entry;
typedef std::list<Vert_decomp_entry>                     Vert_decomp_list;

int main()
{
  // Construct the arrangement.
  Arrangement_2 arr;

  insert_non_intersecting_curve(arr, Segment_2(Point_2(1, 1), Point_2(3, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(1, 1), Point_2(2, 2)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(2, 2), Point_2(3, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(2, 2), Point_2(5, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(3, 2), Point_2(5, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(2, 3), Point_2(3, 3)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(0, 3), Point_2(6, 4)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(4, 4), Point_2(4, 5)));

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
      auto* vh = boost::get<Vertex_const_handle>(&*(curr.first));
      if (vh) std::cout << '(' << (*vh)->point() << ')';
      else {
        auto* hh = boost::get<Halfedge_const_handle>(&*(curr.first));
        if (hh) std::cout << '[' << (*hh)->curve() << ']';
        else {
          auto* fh = boost::get<Face_const_handle>(&*(curr.first));
          CGAL_assertion(fh);
          std::cout << "NONE (" << (*fh)->is_unbounded() << ")";
        }
      }
    }

    std::cout << "   feature above: ";
    if (! curr.second) std::cout << "EMPTY" << std::endl;
    else {
      auto* vh = boost::get<Vertex_const_handle>(&*(curr.second));
      if (vh) std::cout << '(' << (*vh)->point() << ')' << std::endl;
      else {
        auto* hh = boost::get<Halfedge_const_handle>(&*(curr.second));
        if (hh) std::cout << '[' << (*hh)->curve() << ']' << std::endl;
        else {
          auto* fh = boost::get<Face_const_handle>(&*(curr.second));
          CGAL_assertion(fh);
          std::cout << "NONE (" << (*fh)->is_unbounded() << ")" << std::endl;
        }
      }
    }
  }

  return 0;
}


#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <CGAL/graph_traits_dual_arrangement_2.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational                                 Number_type;

typedef CGAL::Cartesian<Number_type>                         Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                   Traits_2;
typedef Traits_2::Point_2                                    Point_2;
typedef Traits_2::X_monotone_curve_2                         Segment_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, unsigned int> Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>                  Arrangement_2;
typedef CGAL::Dual<Arrangement_2>                            Dual_arrangement_2;

int main()
{
  Arrangement_2   arr;

  // Construct an arrangement of seven intersecting line segments.
  insert(arr, Segment_2(Point_2(0, 0), Point_2(0, 6)));
  insert(arr, Segment_2(Point_2(0, 6), Point_2(6, 6)));
  insert(arr, Segment_2(Point_2(6, 6), Point_2(6, 0)));
  insert(arr, Segment_2(Point_2(6, 0), Point_2(0, 0)));

  Arrangement_2::Face_iterator  fit;
  Arrangement_2::Face_handle    f;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (! fit->is_unbounded())
      f = fit;
  }

  insert(arr, Segment_2(Point_2(2, 2), Point_2(2, 4)));
  insert(arr, Segment_2(Point_2(2, 4), Point_2(4, 4)));
  insert(arr, Segment_2(Point_2(4, 4), Point_2(4, 2)));
  insert(arr, Segment_2(Point_2(4, 2), Point_2(2, 2)));
  insert(arr, Segment_2(Point_2(2, 2), Point_2(4, 4)));
  insert(arr, Segment_2(Point_2(4, 2), Point_2(2, 4)));

  Dual_arrangement_2                          darr (arr);
  Dual_arrangement_2::Incident_edge_iterator  fnit = darr.out_edges_begin (f);
  Dual_arrangement_2::Incident_edge_iterator  fnend = darr.out_edges_end (f);
  Dual_arrangement_2::Edge_handle             he;
  unsigned int                                count = 0;

  while (fnit != fnend) {
    he = *fnit;
    ++count;
    std::cout << "Halfedge no. " << count
              << ": (" << he->source()->point()
              << ") --> (" << he->target()->point() << ')' << std::endl;
    ++fnit;
  }
  return 0;
}

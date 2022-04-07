// Testing the handling of isolated vertices.

#include <iostream>

#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;

#define N_POINTS   10

#define TEST_VALIDITY(n)                                                 \
  valid = arr.is_valid();                                                \
  std::cout << "  Inserted " << (n) << " segment(s), arrangement is "    \
            << (valid ? "valid." : "NOT valid!") << std::endl;           \
  if (! valid)                                                           \
  return (1);

int main ()
{
  Arrangement_2   arr;
  Point_2         ps[N_POINTS];
  Vertex_handle   vhs[N_POINTS];
  bool            valid;
  int             k;

  ps[0] = Point_2 (2, 2);
  ps[1] = Point_2 (2, 7);
  ps[2] = Point_2 (4, 9);
  ps[3] = Point_2 (4, 5);
  ps[4] = Point_2 (5, 3);
  ps[5] = Point_2 (7, 1);
  ps[6] = Point_2 (7, 5);
  ps[7] = Point_2 (7, 7);
  ps[8] = Point_2 (9, 3);
  ps[9] = Point_2 (9, 6);

  for (k = 0; k < N_POINTS; k++)
  {
    vhs[k] = insert_point (arr, ps[k]);
  }

  arr.insert_from_left_vertex (Segment_2 (Point_2 (2, 7), Point_2 (4, 7)),
                               vhs[1]);
  TEST_VALIDITY(1);

  arr.insert_from_right_vertex (Segment_2 (Point_2 (6, 6), Point_2 (7, 5)),
                                vhs[6]);
  TEST_VALIDITY(2);

  arr.insert_at_vertices (Segment_2 (Point_2 (7, 1), Point_2 (9, 3)),
                          vhs[5], vhs[8]);
  TEST_VALIDITY(3);

  arr.insert_at_vertices (Segment_2 (Point_2 (7, 5), Point_2 (9, 3)),
                          vhs[6], vhs[8]);
  TEST_VALIDITY(4);

  arr.insert_from_right_vertex (Segment_2 (Point_2 (1, 1), Point_2 (2, 7)),
                                vhs[1]);
  TEST_VALIDITY(5);

  insert_non_intersecting_curve (arr,
                                 Segment_2 (Point_2 (1, 1), Point_2 (7, 1)));
  TEST_VALIDITY(6);

  insert_non_intersecting_curve (arr,
                                 Segment_2 (Point_2 (4, 7), Point_2 (6, 6)));
  TEST_VALIDITY(7);

  insert_non_intersecting_curve (arr,
                                 Segment_2 (Point_2 (2, 7), Point_2 (3, 3)));
  TEST_VALIDITY(8);

  insert_non_intersecting_curve (arr,
                                 Segment_2 (Point_2 (3, 3), Point_2 (7, 1)));
  TEST_VALIDITY(9);

  arr.insert_at_vertices (Segment_2 (Point_2 (7, 5), Point_2 (9, 6)),
                          vhs[6], vhs[9]);
  TEST_VALIDITY(10);

  std::cout << "Arrangement size:"
            << "   V = " << arr.number_of_vertices()
            << " (" << arr.number_of_isolated_vertices() << " isolated)"
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Check the validity more thoroughly.
  valid = is_valid(arr);
  std::cout << "Arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;

  return (0);
}

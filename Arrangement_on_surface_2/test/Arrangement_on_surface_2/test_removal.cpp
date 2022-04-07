// Testing the removal function, covering all possible scenarios.

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
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;

#define N_SEGMENTS 26
#define N_REMOVE   10

Halfedge_handle construct_arr(Arrangement_2& arr)
{
  arr.clear();
  Segment_2 seg;
  seg = Segment_2(Point_2(0, 0), Point_2(2, 0));
  insert_non_intersecting_curve(arr, seg);
  seg = Segment_2(Point_2(2, 4), Point_2(0, 4));
  insert_non_intersecting_curve(arr, seg);
  seg = Segment_2(Point_2(0, 4), Point_2(0, 0));
  insert_non_intersecting_curve(arr, seg);
  seg = Segment_2(Point_2(2, 0), Point_2(4, 0));
  insert_non_intersecting_curve(arr, seg);
  seg = Segment_2(Point_2(4, 0), Point_2(4, 4));
  insert_non_intersecting_curve(arr, seg);
  seg = Segment_2(Point_2(4, 4), Point_2(2, 4));
  insert_non_intersecting_curve(arr, seg);

  seg = Segment_2(Point_2(2, 0), Point_2(2, 4));
  return insert_non_intersecting_curve(arr, seg);
}

/* (1)
 * -------------     -------------
 * |     |     |     |           |
 * |     |  o  |  => |        o  |
 * |     |     |     |           |
 * -------------     -------------
 */
bool test1()
{
  Arrangement_2 arr;
  Halfedge_handle he = construct_arr(arr);
  insert_point(arr, Point_2(3, 2));
  arr.remove_edge(he);
  return arr.is_valid();
}

/* (2)
 * -------------     -------------
 * |  o  |     |     |  o        |
 * |     |  o  |  => |        o  |
 * |  o  |     |     |  o        |
 * -------------     -------------
 */
bool test2()
{
  Arrangement_2 arr;
  Halfedge_handle he = construct_arr(arr);
  insert_point(arr, Point_2(1, 1));
  insert_point(arr, Point_2(1, 3));
  insert_point(arr, Point_2(3, 2));
  arr.remove_edge(he);
  return arr.is_valid();
}

/* (3)
 * -------------     -------------
 * |     |  o  |     |        o  |
 * |  o  |     |  => |  o        |
 * |     |  o  |     |        o  |
 * -------------     -------------
 */
bool test3()
{
  Arrangement_2 arr;
  Halfedge_handle he = construct_arr(arr);
  insert_point(arr, Point_2(1, 2));
  insert_point(arr, Point_2(3, 1));
  insert_point(arr, Point_2(3, 3));
  arr.remove_edge(he);
  return arr.is_valid();
}

/* (4)
 * -------------     -------------
 * |  o  |  o  |     |        o  |
 * |  o  |     |  => |  o        |
 * |  o  |  o  |     |        o  |
 * -------------     -------------
 */
bool test4()
{
  Arrangement_2 arr;
  Halfedge_handle he = construct_arr(arr);
  insert_point(arr, Point_2(1, 1));
  insert_point(arr, Point_2(1, 2));
  insert_point(arr, Point_2(1, 3));
  insert_point(arr, Point_2(3, 1));
  insert_point(arr, Point_2(3, 3));
  arr.remove_edge(he);
  return arr.is_valid();
}

int main ()
{
  // Construct the initial arrangement.
  Arrangement_2   arr;
  Segment_2       segs[N_SEGMENTS];
  Halfedge_handle hhs[N_SEGMENTS];
  bool            valid;
  int             k;

  segs[0] = Segment_2 (Point_2 (5, 9), Point_2 (5, 11));
  segs[1] = Segment_2 (Point_2 (5, 9), Point_2 (7, 9));
  segs[2] = Segment_2 (Point_2 (5, 11), Point_2 (7, 9));
  segs[3] = Segment_2 (Point_2 (7, 6), Point_2 (9, 7));
  segs[4] = Segment_2 (Point_2 (9, 7), Point_2 (7, 9));
  segs[5] = Segment_2 (Point_2 (5, 9), Point_2 (7, 6));
  segs[6] = Segment_2 (Point_2 (10, 11), Point_2 (8, 13));
  segs[7] = Segment_2 (Point_2 (8, 13), Point_2 (11, 12));
  segs[8] = Segment_2 (Point_2 (10, 11), Point_2 (11, 12));
  segs[9] = Segment_2 (Point_2 (1, 20), Point_2 (5, 1));
  segs[10] = Segment_2 (Point_2 (5, 1), Point_2 (12, 6));
  segs[11] = Segment_2 (Point_2 (1, 20), Point_2 (12, 6));
  segs[12] = Segment_2 (Point_2 (13, 13), Point_2 (13, 15));
  segs[13] = Segment_2 (Point_2 (13, 15), Point_2 (15, 12));
  segs[14] = Segment_2 (Point_2 (13, 13), Point_2 (15, 12));
  segs[15] = Segment_2 (Point_2 (11, 12), Point_2 (13, 13));
  segs[16] = Segment_2 (Point_2 (1, 20), Point_2 (11, 17));
  segs[17] = Segment_2 (Point_2 (11, 17), Point_2 (20, 13));
  segs[18] = Segment_2 (Point_2 (12, 6), Point_2 (20, 13));
  segs[19] = Segment_2 (Point_2 (15, 12), Point_2 (20, 13));
  segs[20] = Segment_2 (Point_2 (5, 1), Point_2 (20, 1));
  segs[21] = Segment_2 (Point_2 (20, 1), Point_2 (20, 13));
  segs[22] = Segment_2 (Point_2 (15, 5), Point_2 (17, 3));
  segs[23] = Segment_2 (Point_2 (13, 3), Point_2 (17, 3));
  segs[24] = Segment_2 (Point_2 (12, 6), Point_2 (13, 3));
  segs[25] = Segment_2 (Point_2 (17, 3), Point_2 (20, 1));

  for (k = 0; k < N_SEGMENTS; k++)
  {
    hhs[k] = insert_non_intersecting_curve (arr, segs[k]);
  }
  valid = arr.is_valid();

  std::cout << "Arrangement size:"
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;
  std::cout << "The arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;

  if (! valid)
    return (1);

  // Remove some edges.
  int  del_indices[N_REMOVE] = {25, 23, 22, 1, 3, 11, 19, 15, 4, 24};

  for (k = 0; k < N_REMOVE; k++)
  {
    arr.remove_edge (hhs[del_indices[k]]);
    valid = arr.is_valid();
    std::cout << "  Removed " << k+1 << " segment(s), arrangement is "
              << (valid ? "valid." : "NOT valid!") << std::endl;

    if (! valid)
      return (1);
  }

  std::cout << "Final arrangement size:"
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Check the validity more thoroughly.
  valid = is_valid(arr);
  std::cout << "Arrangement is "
              << (valid ? "valid." : "NOT valid!") << std::endl;

  valid = test1();
  std::cout << "The arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;
  if (!valid) return 1;

  valid = test2();
  std::cout << "The arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;
  if (!valid) return 1;

  valid = test3();
  std::cout << "The arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;
  if (!valid) return 1;

  valid = test4();
  std::cout << "The arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;
  if (!valid) return 1;

  return 0;
}

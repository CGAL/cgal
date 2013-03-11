//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.
// The file ./polyline_example/polylines_illustration.pdf contains
// an illustration of the examples given here.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

#include "arr_print.h"

typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef Segment_traits_2::Curve_2                       Segement_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef Traits_2::X_monotone_curve_2                    X_polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{
  // TODO: Add test with isolated points.
  // TODO: Test with long chains of sub-polylines

  Traits_2 traits;
  Traits_2::Construct_curve_2 polyline_const =
    traits.construct_curve_2_object();
  Traits_2::Construct_x_monotone_curve_2 x_polyline_const =
    traits.construct_x_monotone_curve_2_object();

  // Simple construction test of x-mono polyline from a range of segments.
  std::vector<Segement_2> segs1;
  segs1.push_back(Segement_2(Point_2(0,0),Point_2(1,1)));
  segs1.push_back(Segement_2(Point_2(1,1),Point_2(5,2)));
  X_polyline_2 xpoly1 = x_polyline_const(segs1.begin(),segs1.end());

  Arrangement_2 arr;

  // Polyline from two points.
  // Horizontal line segment.
  Point_2 pt1 = Point_2(-5,2);
  Point_2 pt2 = Point_2(4,2);
  Segement_2 seg = Segement_2(pt1,pt2);
  Polyline_2 poly_one_segment = polyline_const(seg);
  insert(arr,poly_one_segment);

  /* Poyline's outline:
   *
   *     *       *
   *     / \     / \
   *    /   \   /   \
   *   /     \ /     \
   *  *       *       *
   */
  Point_2 points1[5];
  points1[0] = Point_2 (0, 0);
  points1[1] = Point_2 (2, 4);
  points1[2] = Point_2 (3, 0);
  points1[3] = Point_2 (4, 4);
  points1[4] = Point_2 (6, 0);
  Polyline_2 poly1 = polyline_const(&points1[0], &points1[5]);
  insert (arr, poly1);

  // /* Polyline's outline:
  //  *
  //  *      *             *
  //  *     /             / \
  //  *    /             /   \
  //  *   *             *     *
  //  *   \                  /
  //  *    \    *     *     /
  //  *      \ /   \ /   \ /
  //  *       *     *     *
  //  *
  //  */

  // std::list<Point_2>    points2;
  // points2.push_back (Point_2 (1, 3));
  // points2.push_back (Point_2 (0, 2));
  // points2.push_back (Point_2 (1, 0));
  // points2.push_back (Point_2 (2, 1));
  // points2.push_back (Point_2 (3, 0));
  // points2.push_back (Point_2 (4, 1));
  // points2.push_back (Point_2 (5, 0));
  // points2.push_back (Point_2 (6, 2));
  // points2.push_back (Point_2 (5, 3));
  // points2.push_back (Point_2 (4, 2));
  // // Constructing the polyline using the traits class
  // Polyline_2 pi2 = polyline_const(points2.begin(), points2.end());
  // // insert (arr, pi2);

  // /* Polyline's outline:
  //  *
  //  *         *
  //  *        / \
  //  *       /   \
  //  *      /     \
  //  * *---*       *
  //  *
  //  */
  // std::vector<Point_2>  points3 (4);
  // points3[0] = Point_2 (0, 2);
  // points3[1] = Point_2 (1, 2);
  // points3[2] = Point_2 (3, 6);
  // points3[3] = Point_2 (5, 2);
  // Polyline_2 pi3 = polyline_const(points3.begin(), points3.end());
  // // insert (arr, pi3);


  // /*Polyline's outline:
  //  *
  //  *     *
  //  *    / \
  //  *   /   \
  //  *  *     *
  //  *   \   /
  //  *    \ /
  //  *     *
  //  */
  // std::list<Point_2> points4;
  // points4.push_back (Point_2(-2,0));
  // points4.push_back (Point_2(0,-2));
  // points4.push_back (Point_2(2,0));
  // points4.push_back (Point_2(0,2));
  // points4.push_back (Point_2(-2,0));
  // Polyline_2 pi4 = polyline_const(points4.begin(),points4.end());
  // insert (arr, pi4);

  // /* Polyline's outline:
  //  * Axis aligned square
  //  */
  // std::list<Point_2> points5;
  // points5.push_back (Point_2(-2,-2));
  // points5.push_back (Point_2(2,-2));
  // points5.push_back (Point_2(2,2));
  // points5.push_back (Point_2(-2,2));
  // points5.push_back (Point_2(-2,-2));
  // Polyline_2 pi5 = polyline_const(points5.begin(),points5.end());
  // insert (arr, pi5);

  // /* Polyline's outline:
  //  * Single vertical segment
  //  */
  // std::list<Point_2> points6;
  // points6.push_back (Point_2 (0,-7));
  // points6.push_back (Point_2( 0,-4));
  // points6.push_back (Point_2( 0,7));
  // Polyline_2 pi6 = polyline_const(points6.begin(), points6.end());
  // // insert (arr, pi6);

  // /* Polyline's outline:
  //  * Single vertical segment
  //  */
  // std::list<Point_2> points7;
  // points7.push_back (Point_2 (-4,-4));
  // points7.push_back (Point_2 (-2,-4));
  // points7.push_back (Point_2 (-1,-4));
  // points7.push_back (Point_2 (-1,2));
  // points7.push_back (Point_2 (-1,3));
  // points7.push_back (Point_2 (1,3));
  // points7.push_back (Point_2 (2,3));
  // points7.push_back (Point_2 (2,4));
  // Polyline_2 pi7 = polyline_const(points7.begin(), points7.end());
  // // insert (arr, pi7);

  // Print the arrangement
  print_arrangement (arr);
  return 0;
}

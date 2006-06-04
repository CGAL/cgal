//! \file examples/Arrangement_2/ex_infinite_non_intersecting.C
// Constructing an arrangement of unbounded linear objects using the insertion
// function for non-intersecting curves.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

typedef int                                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Segment_2                           Segment_2;
typedef Traits_2::Ray_2                               Ray_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;

int main ()
{
  // Construct 3 rays and 6 line segments.
  const int          n_curves = 9;
  X_monotone_curve_2 curves[9];

  curves[0] = Ray_2 (Point_2 (-5, -1), Point_2 (-6, -2));
  curves[1] = Segment_2 (Point_2 (-5, -1), Point_2 (5, -1));
  curves[2] = Segment_2 (Point_2 (-5, -1), Point_2 (0, 8));
  curves[3] = Segment_2 (Point_2 (-1, 2), Point_2 (1, 2));
  curves[4] = Segment_2 (Point_2 (-1, 2), Point_2 (0, 0));
  curves[5] = Segment_2 (Point_2 (0, 0), Point_2 (1, 2));
  curves[6] = Ray_2 (Point_2 (0, 8), Point_2 (0, 9));
  curves[7] = Segment_2 (Point_2 (0, 8), Point_2 (5, -1));
  curves[8] = Ray_2 (Point_2 (5, -1), Point_2 (6, -2));

  // Construct the arrangement by inserting the curves incermentally.
  Arrangement_2      arr;
  Naive_pl           naive_pl (arr);
  int                k;

  for (k = 0; k < n_curves; k++)
    insert_non_intersecting_curve (arr, curves[k], naive_pl);

  // Print out the size of the resulting arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (" << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return (0);
}

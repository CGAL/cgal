//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <list>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // Construct the arrangement of five intersecting segments.
  Arrangement_2           arr;
  std::list<Segment_2>    segments;

  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(2, 4)));
  segments.push_back (Segment_2 (Point_2(5, 0), Point_2(5, 5)));
  segments.push_back (Segment_2 (Point_2(1, 0), Point_2(5, 3)));  
  segments.push_back (Segment_2 (Point_2(0, 2), Point_2(6, 0)));
  segments.push_back (Segment_2 (Point_2(3, 0), Point_2(5, 5)));

  insert (arr, segments.begin(), segments.end());

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}

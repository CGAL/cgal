#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               Coord_type;
typedef K::Point_2                                          Point;
typedef CGAL::Delaunay_triangulation_2<K>                   Delaunay_triangulation;

// Resulting points-coordinates pairs will be stored in an object of this type
typedef std::vector<std::pair<Point, Coord_type> >          Point_coordinate_vector;

int main()
{
  Delaunay_triangulation dt;

  for(int y=0; y<3; ++y)
    for(int x=0; x<3; ++x)
      dt.insert(K::Point_2(x, y));

  // coordinates computation
  K::Point_2 p(1.2, 0.7); // query point
  Point_coordinate_vector coords;

  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, K::FT, bool> result =
      CGAL::natural_neighbor_coordinates_2(dt, p, std::back_inserter(coords));

  if(!result.third)
  {
    std::cout << "The coordinate computation was not successful." << std::endl;
    std::cout << "The point (" << p << ") lies outside the convex hull." << std::endl;
  }

  K::FT norm = result.second;
  std::cout << "Coordinate computation successful." << std::endl;
  std::cout << "Normalization factor: " << norm << std::endl;

  std::cout << "Coordinates for point: (" << p << ") are the following: " << std::endl;
  for(std::size_t i=0; i<coords.size(); ++i)
    std::cout << "  Point: (" << coords[i].first << ") coeff: " << coords[i].second << std::endl;

  return EXIT_SUCCESS;
}

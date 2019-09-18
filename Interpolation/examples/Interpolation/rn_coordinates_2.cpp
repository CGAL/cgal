#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               FT;
typedef CGAL::Regular_triangulation_2<K>                    Regular_triangulation;
typedef Regular_triangulation::Bare_point                   Bare_point;
typedef Regular_triangulation::Weighted_point               Weighted_point;

typedef std::vector<std::pair<Weighted_point, FT> >         Point_coordinate_vector;

int main()
{
  Regular_triangulation rt;

  for (int y=0; y<3; ++y)
    for (int x=0; x<3; ++x)
      rt.insert(Weighted_point(Bare_point(x, y), 0. /*weight*/));

  // coordinate computation
  Weighted_point wp(Bare_point(1.2, 0.7), 2.);
  Point_coordinate_vector coords;
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, K::FT, bool> result =
      CGAL::regular_neighbor_coordinates_2(rt, wp, std::back_inserter(coords));

  if(!result.third)
  {
    std::cout << "The coordinate computation was not successful." << std::endl;
    std::cout << "The point (" <<wp.point() << ") lies outside the convex hull." << std::endl;
  }

  K::FT norm = result.second;
  std::cout << "Coordinate computation successful." << std::endl;
  std::cout << "Normalization factor: " << norm << std::endl;

  std::cout << "Coordinates for point: (" << wp << ") are the following: " << std::endl;
  for(std::size_t i=0; i<coords.size(); ++i)
    std::cout << "  Point: (" << coords[i].first << ") coeff: " << coords[i].second << std::endl;

  return EXIT_SUCCESS;
}

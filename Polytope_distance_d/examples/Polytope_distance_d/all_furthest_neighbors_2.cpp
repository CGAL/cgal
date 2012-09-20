#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <iostream>
#include <vector>

typedef double                                    FT;

typedef CGAL::Cartesian<FT>                       Kernel;

typedef Kernel::Point_2                           Point;
typedef std::vector<int>                          Index_cont;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point>    Generator;
typedef CGAL::Ostream_iterator<int,std::ostream>  Oiterator;

int main()
{
  // generate random convex polygon:
  Polygon_2 p;
  CGAL::random_convex_set_2(10, std::back_inserter(p), Generator(1));

  // compute all furthest neighbors:
  CGAL::all_furthest_neighbors_2(p.vertices_begin(), p.vertices_end(),
                                 Oiterator(std::cout));
  std::cout << std::endl;

  return 0;
}

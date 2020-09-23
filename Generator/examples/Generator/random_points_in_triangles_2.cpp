#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>

#include <iostream>
#include <fstream>
using namespace CGAL;
typedef Simple_cartesian<double>                           K;
typedef K::Point_2                                         Point;


int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  // Create input triangles
  std::vector<K::Triangle_2> triangles;
  for(int i=0; i< 5; ++i)
  {
    triangles.push_back(K::Triangle_2(Point(i,0), Point(i+1,0), Point(i+0.5,1)));
  }

  // Create the generator, input is the vector of Triangle_2
  Random_points_in_triangles_2<Point> g(triangles);
  // Get 100 random points in cdt
  std::copy_n(g, 1000, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 1000);

  // print the first point that was generated
  std::cout << points[0] << std::endl;

  return 0;
}



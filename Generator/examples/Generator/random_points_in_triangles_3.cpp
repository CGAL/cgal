#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>
using namespace CGAL;
typedef Simple_cartesian<double>                           K;
typedef K::Point_3                                         Point;


int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  // Create input triangles
  std::vector<K::Triangle_3> triangles;
  for(int i=0; i< 5; ++i)
  {
    triangles.push_back(K::Triangle_3(Point(i,0,0.5*i), Point(i+1,0,0.5*i), Point(i+0.5,1,0.5*i)));
  }

  // Create the generator, input is the vector of Triangle_3
  Random_points_in_triangles_3<Point> g(triangles);
  // Get 100 random points in cdt
  std::copy_n(g, 1000, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 1000);

  // print the first point that was generated
  std::cout << points[0] << std::endl;

  return 0;
}



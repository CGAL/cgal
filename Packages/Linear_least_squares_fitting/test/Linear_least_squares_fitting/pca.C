// Test program for the linear_least_square_fitting() functions.
// Pierre Alliez

#include <vector>
#include <cassert>
#include <stdlib.h>

#include <CGAL/Cartesian.h>

#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Point_2           Point_2;
typedef K::Line_2           Line_2;
typedef K::Point_3           Point_3;

Point_2 random_point_2()
{
  FT x = rand() / (FT)RAND_MAX;
  FT y = rand() / (FT)RAND_MAX;
  return Point_2(x,y);
}

void test_2(const unsigned int nb_points)
{
  std::cout << "2D case" << std::endl;

  // create points that are collinear to two 
  // randomly chosen points.
  std::vector<Point_2> points;
  points.push_back(random_point_2());
  points.push_back(random_point_2());
  std::cout << "  generate " << nb_points << 
       " random points on a segment...";
  random_collinear_points_2(points.begin(),points.end(),nb_points,
			    std::back_inserter(points));
  std::cout << "done" << std::endl;

  // fit a line
  std::cout << "  fit a line...";
  Line_2 line;
  linear_least_squares_fitting_2(points.begin(),points.end(),line);
  std::cout << "done" << std::endl;
  // todo: output in a ps file, add asserts, etc
}

int main()
{
  test_2(100);
  return 0;
}

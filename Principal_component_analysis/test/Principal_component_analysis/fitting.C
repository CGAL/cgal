// Demo program for the linear_least_square_fitting() functions.
#include <vector>
#include <cassert>
#include <stdlib.h>

#include <CGAL/Cartesian.h>

#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

// 2D types
typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Triangle_2 Triangle_2;

// 3D types
typedef Kernel::Line_3 Line_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Triangle_3 Triangle_3;

void test_2_point_set(const unsigned int nb_points)
{
  std::cout << "2D: fit a line to a point set" << std::endl;

  // create random points on a horizontal segment
  std::vector<Point_2> points;
  Point_2 p = Point(0.0,0.0);
  Point_2 q = Point(1.0,0.0);
  std::cout << "  generate " << nb_points << 
       " 2D random points on a unit horizontal segment...";
  points_on_segment_2(p,q,nb_points,std::back_inserter(points));
  std::cout << "done" << std::endl;

  // fit a line
  std::cout << "  fit a 2D line...";
  Line_2 line;
  FT quality = linear_least_squares_fitting_2(points.begin(),points.end(),line);
  std::cout << "done (quality: " << quality << ")" << std::endl;

  // TODO: check fitting line is ~horizontal
  //if(!horizontal)
  //exit(1); // failure
}


int main()
{
  std::cout << "Test linear_least_squares_fitting"  << std::endl;

  test_2D_point_set(100);
  test_3D_point_set(100);
  // test_3D_triangle_set(100);

  return 0; // success
}

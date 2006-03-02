// test for the linear_least_square_fitting() functions.
#include <vector>
#include <cassert>
#include <stdlib.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_3.h>

// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Triangle_3 Triangle;

// case with only one point in container
// the fitting plane must be horizontal by default
void test_3D()
{
  std::list<Point_3> points;
  points.push_back(Point_3(0,0,0));

  // fit a plane
  // call all versions of the function
  std::cout << "fit 3D plane...";
  Kernel k;
  Plane plane;
  Point centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane);
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,centroid);
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,centroid,k);
  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!plane.is_horizontal())
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}


// case with a point set on a horizontal plane
// the fitting plane must be horizontal
void test_3D_point_set(const unsigned int nb_points)
{
}


int main()
{
  std::cout << "Test linear_least_squares_fitting"  << std::endl;

  // 3D
  test_3D();
  test_3D_point_set(100);
  // test_3D_triangle_set(100);

  return 0; // success
}

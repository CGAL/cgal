// test for the linear_least_square_fitting() functions.
#include <vector>
#include <cassert>
#include <stdlib.h>

#include <CGAL/Cartesian.h>
#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/point_generators_2.h>

// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Vector_2 Vector_2;

typedef Kernel::Line_3 Line_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Triangle_3 Triangle_3;

// case with only one point in container
// the fitting line must be horizontal by default
void test_2D()
{
  std::vector<Point_2> points;
  points.push_back(Point_2(0.0,0.0));

  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line...";
  Kernel k;
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line);
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid);
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,k);
  std::cout << "done (quality: " << quality << ")" << std::endl;

  Vector_2 v = line.to_vector();
  v = v / std::sqrt(v*v); // normalize
  std::cout << "vector: " << v.x() << " " 
                          << v.y() << std::endl;

  if(!line.is_horizontal())
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}


// case with a point set on a horizontal segment
// the fitting line must be horizontal
void test_2D_point_set(const unsigned int nb_points)
{
  // create random points on a horizontal segment
  std::vector<Point_2> points;
  points.push_back(Point_2(0.0,0.0));
  points.push_back(Point_2(1.0,0.0));
  std::cout << "generate " << nb_points << 
       " 2D random collinear points on a horizontal line...";
  random_collinear_points_2(points.begin(),points.end(),nb_points,
                            std::back_inserter(points));
  std::cout << "done " << std::endl;

  // fit a line
  std::cout << "fit 2D line...";
  Line_2 line;
  Point_2 centroid;

  // call all versions of the function
  FT quality;
  Kernel k;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line);
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid);
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,k);

  std::cout << "done (quality: " << quality << ")" << std::endl;

  FT y1 = line.y_at_x((FT)0);
  FT y2 = line.y_at_x((FT)1);

  std::cout << "y1" << y1 << " " << y2 << std::endl;

  if(!line.is_horizontal())
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}


int main()
{
  std::cout << "Test linear_least_squares_fitting"  << std::endl;

  // 2D
  test_2D();
  test_2D_point_set(100);

  // 3D
  // test_3D();
  // test_3D_point_set(100);
  // test_3D_triangle_set(100);

  return 0; // success
}

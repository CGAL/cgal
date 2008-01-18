// test for the linear_least_square_fitting() functions.


#include <CGAL/Cartesian.h>
#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>

#include <vector>
#include <CGAL/Testsuite/assert.h>
#include <cstdlib>

// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Vector_2 Vector_2;

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
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::PCA_dimension_0_tag(),k);
  std::cout << "done (quality: " << quality << ")" << std::endl;

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
  // create points on a horizontal segment
  Point_2 p(0.0,0.5);
  Point_2 q(1.0,0.5);

  std::cout << "generate " << nb_points << 
       " 2D points on a horizontal line...";
  std::list<Point_2> points;
  points_on_segment_2(p,q,100,std::back_inserter(points));
  std::cout << "done " << std::endl;

  // fit a line
  std::cout << "fit 2D line...";
  Line_2 line;
  Point_2 centroid;

  // call all versions of the function
  FT quality;
  Kernel k;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::PCA_dimension_0_tag(),k);

  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!line.is_horizontal())
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}


int main()
{
  std::cout << "Test 2D linear_least_squares_fitting"  << std::endl;

  test_2D();
  test_2D_point_set(100);

  return 0; // success
}

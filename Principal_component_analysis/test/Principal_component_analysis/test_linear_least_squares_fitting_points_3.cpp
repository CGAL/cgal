#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

// test for the linear_least_square_fitting() functions.


#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Default_diagonalize_traits.h>

#include <cassert>

#include <cstdlib>
#include <list>

// types
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;

void test_one_point();
void test_point_set(const unsigned int nb_points);


int main(void)
{
  std::cout << "Test linear least squares fitting for 3D points"  << std::endl;
  test_one_point();
  test_point_set(1000);
  return 0; // success
}


void fit_point_set(std::list<Point>& points,
                   Plane& plane,
                   Line& line)
{
  // fit a plane
  // call all versions of the function
  Kernel kernel;
  FT quality;
  Point centroid;

  std::cout << "fit 3D line...";
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),line,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>(),kernel,
                                           CGAL::Default_diagonalize_traits<FT,3>());

  std::cout << "done (quality: " << quality << ")" << std::endl;

  std::cout << "fit 3D plane...";
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,centroid,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,centroid,CGAL::Dimension_tag<0>(),kernel,
                                           CGAL::Default_diagonalize_traits<FT,3>());

  std::cout << "done (quality: " << quality << ")" << std::endl;
}

// case with only one point in container
// the fitting plane must be horizontal by default
void test_one_point()
{
  std::list<Point> points;
  points.push_back(Point(FT(0),FT(0),FT(0)));

  // fit plane
  Plane plane;
  Line line;
  fit_point_set(points,plane,line);

  Point point(FT(0),FT(0),FT(0));
  Vector vec(FT(0),FT(0),FT(1));
  Plane horizontal_plane(point,vec);
  if(!parallel(horizontal_plane,plane))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}

Point random_point_xy()
{
  FT x = (FT)((double)std::rand() / (double)RAND_MAX);
  FT y = (FT)((double)std::rand() / (double)RAND_MAX);
  return Point(x,y,FT(0));
}

// case with a random point set on a horizontal plane
// the fitting plane must be horizontal
void test_point_set(const unsigned int nb_points)
{
  std::list<Point> points;
  unsigned int i;
  for(i=0;i<nb_points;i++)
    points.push_back(random_point_xy());

  // fit plane
  Plane plane;
  Line line;
  fit_point_set(points,plane,line);

  Point point(FT(0),FT(0),FT(0));
  Vector vec(FT(0),FT(0),FT(1));
  Plane horizontal_plane(point,vec);
  if(!parallel(horizontal_plane,plane))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}


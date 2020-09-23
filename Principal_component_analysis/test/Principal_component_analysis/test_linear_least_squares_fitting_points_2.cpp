#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

// test for the linear_least_square_fitting() functions.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/algorithm.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Default_diagonalize_traits.h>

#include <vector>
#include <cassert>
#include <cstdlib>

// types

// case with only one point in container
// the fitting line must be horizontal by default
template <typename Kernel>
void test_2D()
{
        typedef typename Kernel::FT FT;
        typedef typename Kernel::Line_2 Line_2;
        typedef typename Kernel::Point_2 Point_2;

  std::vector<Point_2> points;
  points.push_back(Point_2(FT(0),FT(0)));

  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line...";
  Kernel k;
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>(),k,
                                           CGAL::Default_diagonalize_traits<FT,2>());

  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!line.is_horizontal())
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}


// case with a point set on a horizontal segment
// the fitting line must be horizontal
template <typename Kernel>
void test_2D_point_set(const unsigned int nb_points)
{
        typedef typename Kernel::FT FT;
        typedef typename Kernel::Line_2 Line_2;
        typedef typename Kernel::Point_2 Point_2;

  // create points on a horizontal segment
  Point_2 p(FT(0.0),FT(0.5));
  Point_2 q(FT(1.0),FT(0.5));

  std::cout << "generate " << nb_points <<
       " 2D points on a horizontal line...";
  std::list<Point_2> points;
  points_on_segment_2(p, q, 100, std::back_inserter(points));
  std::cout << "done " << std::endl;

  // fit a line
  std::cout << "fit 2D line...";
  Line_2 line;
  Point_2 centroid;

  // call all versions of the function
  Kernel k;
  FT quality;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>(),k,
                                           CGAL::Default_diagonalize_traits<FT,2>());

  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!line.is_horizontal())
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}


int main()
{
  std::cout << "Test 2D linear least squares fitting of points"  << std::endl;

        typedef CGAL::Simple_cartesian<double> Kernel_double;
  test_2D<Kernel_double>();
  test_2D_point_set<Kernel_double>(100);

        typedef CGAL::Simple_cartesian<float> Kernel_float;
  test_2D<Kernel_float>();
  test_2D_point_set<Kernel_float>(100);

  return 0; // success
}

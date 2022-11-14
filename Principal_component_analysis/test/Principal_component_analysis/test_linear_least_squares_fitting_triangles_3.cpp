#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

// Example program for the linear_least_square_fitting function on set of triangles in 3D
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Diagonalize_traits.h>
#endif

#include <list>

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  Kernel;
typedef Kernel::Line_3              Line;
typedef Kernel::Plane_3             Plane;
typedef Kernel::Point_3             Point;
typedef Kernel::Triangle_3          Triangle;

int main()
{
  std::cout << "Test 3D linear least squares fitting of triangles"  << std::endl;

  std::list<Triangle> triangles;
        Point a(0.0,0.0,0.0);
        Point b(1.0,0.0,0.0);
        Point c(0.0,1.0,0.0);
        Point d(0.0,0.0,1.0);
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));

  Line line;
  Plane plane;
  Point centroid;

        // no centroid, fit line and plane
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<0>());

        // centroid, fit line and plane
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid, CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::Dimension_tag<0>());

  // If Eigen is available, it's used by default everywhere.
  // These additional lines test the fallback version
#ifdef CGAL_EIGEN3_ENABLED
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid, CGAL::Dimension_tag<2>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid, CGAL::Dimension_tag<1>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid, CGAL::Dimension_tag<0>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::Dimension_tag<2>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::Dimension_tag<1>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::Dimension_tag<0>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
#endif

  return 0;
}

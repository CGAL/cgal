#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

// Example program for the linear_least_square_fitting function on set of cuboids in 3D
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Diagonalize_traits.h>
#endif

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  Kernel;
typedef Kernel::Line_3              Line;
typedef Kernel::Plane_3             Plane;
typedef Kernel::Point_3             Point;
typedef Kernel::Iso_cuboid_3        Iso_cuboid;

int main()
{
  std::cout << "Test linear least squares fitting of 3D cuboids"  << std::endl;
  std::list<Iso_cuboid> cuboids;
        Point a(0.0,0.0,0.0);
        Point b(1.0,2.0,3.0);
        Point c(4.0,5.0,6.0);
  cuboids.push_back(Iso_cuboid(a,b));
  cuboids.push_back(Iso_cuboid(a,c));

  Line line;
  Plane plane;
  Point centroid;

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::Dimension_tag<0>());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<0>());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<0>());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<0>());

  // If Eigen is available, it's used by default everywhere.
  // These additional lines test the fallback version
#ifdef CGAL_EIGEN3_ENABLED
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<3>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<2>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<1>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::Dimension_tag<0>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<3>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<2>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<1>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::Dimension_tag<0>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
#endif

  return 0;
}

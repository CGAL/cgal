// Example program for the linear_least_square_fitting function

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_diagonalize_traits.h>
#endif
#include <CGAL/Diagonalize_traits.h>

#include <list>

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  Kernel;
typedef Kernel::Line_3              Line;
typedef Kernel::Plane_3             Plane;
typedef Kernel::Point_3             Point;
typedef Kernel::Sphere_3            Sphere;

int main(void)
{
  std::cout << "Test linear least squares fitting of 3D spheres"  << std::endl;

        // centers
        Point c1(0.0,0.0,0.0);
        Point c2(1.0,1.0,1.0);

        // radii
        FT sqr1 = 0.1;
        FT sqr2 = 0.5;

        // add two spheres
  std::list<Sphere> spheres;
  spheres.push_back(Sphere(c1,sqr1));
  spheres.push_back(Sphere(c2,sqr2));

  Line line;
  Plane plane;
  Kernel kernel;
  Point centroid;


  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,centroid,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,centroid,CGAL::Dimension_tag<3>(),kernel,
                                 CGAL::Diagonalize_traits<FT,3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,centroid,CGAL::Dimension_tag<2>(),kernel,
                                 CGAL::Diagonalize_traits<FT,3>());
#ifdef CGAL_EIGEN3_ENABLED
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,centroid,CGAL::Dimension_tag<3>(),kernel,
                                 CGAL::Eigen_diagonalize_traits<FT,3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,centroid,CGAL::Dimension_tag<2>(),kernel,
                                 CGAL::Eigen_diagonalize_traits<FT,3>());
#endif
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,centroid,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,centroid,CGAL::Dimension_tag<3>(),kernel,
                                 CGAL::Diagonalize_traits<FT,3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,centroid,CGAL::Dimension_tag<2>(),kernel,
                                 CGAL::Diagonalize_traits<FT,3>());
#ifdef CGAL_EIGEN3_ENABLED
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,centroid,CGAL::Dimension_tag<3>(),kernel,
                                 CGAL::Eigen_diagonalize_traits<FT,3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,centroid,CGAL::Dimension_tag<2>(),kernel,
                                 CGAL::Eigen_diagonalize_traits<FT,3>());
#endif

  return 0;
}

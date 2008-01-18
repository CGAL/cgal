// Example program for the linear_least_square_fitting function on set of cuboids in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  Kernel;
typedef Kernel::Line_3       Line;
typedef Kernel::Plane_3      Plane;
typedef Kernel::Point_3      Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;

int main()
{
  std::cout << "Test 3D linear least squares fitting of cuboids"  << std::endl;

  std::list<Iso_cuboid> cuboids;
	Point a(0.0,0.0,0.0);
	Point b(1.0,2.0,3.0);
	Point c(4.0,5.0,6.0);
  cuboids.push_back(Iso_cuboid(a,b));
  cuboids.push_back(Iso_cuboid(a,c));

  Kernel kernel;
  Line line;
  Plane plane;
  Point centroid;

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_0_tag());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,centroid,CGAL::PCA_dimension_0_tag());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_0_tag());

  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,centroid,CGAL::PCA_dimension_0_tag());

  return 0;
}

// Example program for the linear_least_square_fitting function

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  Kernel;
typedef Kernel::Line_3       Line_3;
typedef Kernel::Plane_3      Plane_3;
typedef Kernel::Point_3      Point_3;
typedef Kernel::Sphere_3     Sphere_3;

int main(void)
{
	// centers
	Point c1(0.0,0.0,0.0);
	Point c2(1.0,1.0,1.0);

	// radii
	FT sqr1 = 0.1;
	FT sqr2 = 0.5;

	// add two spheres
  std::list<Sphere_3> spheres;
  spheres.push_back(Sphere_3(c1,sqr1));
  spheres.push_back(Sphere_3(c2,sqr2));

  Kernel k;
  Point_3 c;
  Line_3 line;
  Plane_3 plane;

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,c,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,c,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,k,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,k,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,c,k,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,c,k,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,c,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,c,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,k,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,k,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,c,k,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,c,k,CGAL::PCA_dimension_2_tag());


  return 0;
}

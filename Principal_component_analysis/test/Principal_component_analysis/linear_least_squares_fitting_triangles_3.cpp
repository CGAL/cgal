// Example program for the linear_least_square_fitting function on set of triangles in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  Kernel;
typedef Kernel::Line_3       Line;
typedef Kernel::Plane_3      Plane;
typedef Kernel::Point_3      Point;
typedef Kernel::Triangle_3   Triangle;

int main()
{
  std::list<Triangle> triangles;

	Point a(0.0,0.0,0.0);
	Point b(1.0,0.0,0.0);
	Point c(0.0,1.0,0.0);
	Point d(0.0,0.0,1.0);
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));

  Line line;
  Plane plane;
	Kernel kernel;
  Point centroid;

  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_0_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid,CGAL::PCA_dimension_0_tag());

  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_0_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::PCA_dimension_0_tag());
  
  return 0;
}

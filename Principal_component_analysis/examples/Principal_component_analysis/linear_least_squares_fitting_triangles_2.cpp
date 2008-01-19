// Example program for linear least squares fitting of 2D triangles
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line;
typedef K::Point_2           Point;
typedef K::Triangle_2        Triangle;

int main()
{
	// generate 2D triangles
  std::list<Triangle> triangles;
	Point a(1.0,2.0,3.0);
	Point b(4.0,5.0,6.0);
	Point c(7.0,8.0,9.0);
	Point d(0.1,0.2,0.3);
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));

  // fit line to triangles
  Line line;
  linear_least_squares_fitting_2(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_2_tag());

  // fit line to triangle edges
  linear_least_squares_fitting_2(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_1_tag());

	// fit line to triangle vertices
  linear_least_squares_fitting_2(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_0_tag());

  return 0;
}

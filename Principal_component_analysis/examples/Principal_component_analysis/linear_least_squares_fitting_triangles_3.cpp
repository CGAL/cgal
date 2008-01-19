// Example program for the linear_least_square_fitting function on set of triangles in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Triangle_3        Triangle;

int main(void)
{
  std::list<Triangle> triangles;
	Point a( 0.0,0.0,0.0);
	Point b( 1.0,0.0,0.0);
	Point c(-1.0,0.0,0.0);
	Point d( 0.0,1.0,1.0);
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));

  Line line;
  Plane plane;

	// fit a line and a plane to triangles
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_2_tag());
  
	// fit a line and a plane to triangle edges
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_2_tag());
  
	// fit a line and a plane to triangle vertices
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::PCA_dimension_0_tag());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_0_tag());
  
  return 0;
}

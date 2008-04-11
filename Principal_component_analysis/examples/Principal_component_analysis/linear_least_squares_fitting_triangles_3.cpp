// Example program for the linear_least_square_fitting function on set of triangles in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>
#include <cstdlib> // for std::rand

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Triangle_3        Triangle;

FT random_value()
{
	return (FT)std::rand() / (FT)RAND_MAX;
}

Point random_point()
{
	return Point(random_value(),
		           random_value(),
							 random_value());
}

int main(void)
{
  std::list<Triangle> triangles;
	Point a = random_point();
	Point b = random_point();
	Point c = random_point();
	Point d = random_point();
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));

  Line line;
  Plane plane;

	// fit a line and a plane to triangles
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
  
	// fit a line and a plane to triangle edges
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
  
	// fit a line and a plane to triangle vertices
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<0>());
  
  return 0;
}

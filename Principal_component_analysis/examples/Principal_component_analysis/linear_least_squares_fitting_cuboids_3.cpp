// Example program for linear least squares fitting of 3D cuboids
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Iso_cuboid_3      Iso_cuboid;

int main()
{
  std::list<Iso_cuboid> cuboids;
	Point a(1.0,2.0,3.0);
	Point b(4.0,5.0,6.0);
	Point c(7.0,8.0,9.0);
  cuboids.push_back(Iso_cuboid(a,b));
  cuboids.push_back(Iso_cuboid(a,c));

  Line line;
  Plane plane;

	// fit volumes
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line, CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<3>());
  
  // fit quadrangle faces
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<2>());

	// fit edges
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line, CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<1>());

	// fit vertices
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::Dimension_tag<0>());

  return 0;
}

// Example program for linear least squares fitting of 3D spheres and balls
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Sphere_3          Sphere;

int main()
{
  std::list<Sphere> spheres;
  spheres.push_back(Sphere(Point(1.0,2.0,3.0),16.0));
  spheres.push_back(Sphere(Point(4.0,5.0,6.0),25.0));

  Line line;
  Plane plane;

  // fit a line and a plane to spheres (dimension 2)
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::Dimension_tag<2>());

  // fit a line and a plane to balls (dimension 3)
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line, CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::Dimension_tag<3>());
  
  return 0;
}

// Example program for the linear_least_square_fitting function on set of triangles in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Plane_3           Plane_3;
typedef K::Point_3           Point_3;
typedef K::Triangle_3        Triangle_3;

int main(void)
{
  std::list<Triangle_3> triangles;
  triangles.push_back(Triangle_3(Point_3(1.0,0.0,0.0),Point_3(0.0,1.0,0.0),Point_3(0.0,0.0,0.0)));
  triangles.push_back(Triangle_3(Point_3(-1.0,0.0,0.0),Point_3(0.0,-1.0,0.0),Point_3(0.0,0.0,0.0)));

  Line_3 line;
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_2_tag());

  Plane_3 plane;
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_2_tag());
  
  return 0;
}

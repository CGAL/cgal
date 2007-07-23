// Example program for the linear_least_square_fitting function on set of tetrahedrons in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Plane_3           Plane_3;
typedef K::Point_3           Point_3;
typedef K::Tetrahedron_3          Tetrahedron_3;

int main()
{
  std::list<Tetrahedron_3> tetrahedrons;
  tetrahedrons.push_back(Tetrahedron_3(Point_3(0.0,0.0,0.0),Point_3(1.0,0.0,0.0),Point_3(0.0,1.0,0.0),Point_3(0.0,0.0,1.0)));
  tetrahedrons.push_back(Tetrahedron_3(Point_3(0.0,0.0,0.0),Point_3(-1.0,0.0,0.0),Point_3(0.0,-1.0,0.0),Point_3(0.0,0.0,1.0)));

  Line_3 line;
  linear_least_squares_fitting_3(tetrahedrons.begin(),tetrahedrons.end(),line,CGAL::PCA_dimension_3_tag());

  Plane_3 plane;
  linear_least_squares_fitting_3(tetrahedrons.begin(),tetrahedrons.end(),plane,CGAL::PCA_dimension_3_tag());
    
  return 0;
}

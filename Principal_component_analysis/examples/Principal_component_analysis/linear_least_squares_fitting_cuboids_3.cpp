// Example program for the linear_least_square_fitting function on set of cuboids in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Plane_3          Plane_3;
typedef K::Point_3           Point_3;
typedef K::Iso_cuboid_3   Iso_cuboid_3;

int main()
{
  std::list<Iso_cuboid_3> cuboids;
  cuboids.push_back(Iso_cuboid_3(Point_3(0.0,0.0,0.0),Point_3(1.0,1.0,1.0)));
  cuboids.push_back(Iso_cuboid_3(Point_3(1.0,1.0,1.0),Point_3(6.0,6.0,6.0)));

  Line_3 line;
  Plane_3 plane;

	// fit volume
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_3_tag());
  
  // fit faces
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_2_tag());

	// fit edges
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_1_tag());

	// fit vertices
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),line,CGAL::PCA_dimension_0_tag());
  linear_least_squares_fitting_3(cuboids.begin(),cuboids.end(),plane,CGAL::PCA_dimension_0_tag());

  return 0;
}

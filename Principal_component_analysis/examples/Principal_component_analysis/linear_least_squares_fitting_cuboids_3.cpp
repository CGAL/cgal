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
  std::list<Iso_cuboid_3> Iso_cuboids;
  Iso_cuboids.push_back(Iso_cuboid_3(Point_3(0.0,0.0,0.0),Point_3(1.0,1.0,1.0)));
  Iso_cuboids.push_back(Iso_cuboid_3(Point_3(1.0,1.0,1.0),Point_3(6.0,6.0,6.0)));

  Line_3 line;
  linear_least_squares_fitting_3(Iso_cuboids.begin(),Iso_cuboids.end(),line,CGAL::PCA_dimension_3_tag());

  Plane_3 plane;
  linear_least_squares_fitting_3(Iso_cuboids.begin(),Iso_cuboids.end(),plane,CGAL::PCA_dimension_3_tag());
  

  //Use only faces.
  linear_least_squares_fitting_3(Iso_cuboids.begin(),Iso_cuboids.end(),line,CGAL::PCA_dimension_2_tag());

  linear_least_squares_fitting_3(Iso_cuboids.begin(),Iso_cuboids.end(),plane,CGAL::PCA_dimension_2_tag());

  return 0;
}

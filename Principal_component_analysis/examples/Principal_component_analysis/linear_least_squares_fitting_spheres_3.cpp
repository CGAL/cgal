// Example program for the linear_least_square_fitting function

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Plane_3           Plane_3;
typedef K::Point_3           Point_3;
typedef K::Sphere_3          Sphere_3;

int main()
{
  std::list<Sphere_3> spheres;
  spheres.push_back(Sphere_3(Point_3(0.0,0.0,0.0),9));
  spheres.push_back(Sphere_3(Point_3(0.0,10.0,0.0),25));

  Line_3 line;
  FT i = linear_least_squares_fitting_3(spheres.begin(),spheres.end(),line,CGAL::PCA_dimension_3_tag());

  std::cout<<"Line's accuracy: "<<i<<std::endl;
  std::cout<<"Line: "<<line<<std::endl;

  Plane_3 plane;
  FT j = linear_least_squares_fitting_3(spheres.begin(),spheres.end(),plane,CGAL::PCA_dimension_3_tag());

  std::cout<<"Plane's accuracy: "<<j<<std::endl;
  std::cout<<"Plane: "<<plane<<std::endl;
  return 0;
}

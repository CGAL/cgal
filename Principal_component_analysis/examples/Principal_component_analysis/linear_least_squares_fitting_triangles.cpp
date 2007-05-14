// Example program for the linear_least_square_fitting function

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line_2;
typedef K::Point_2           Point_2;
typedef K::Triangle_2        Triangle_2;

int main()
{
  std::list<Triangle_2> triangles;
  triangles.push_back(Triangle_2(Point_2(1.0,0.0),Point_2(0.0,1.0),Point_2(1.0,1.0)));
  triangles.push_back(Triangle_2(Point_2(-1.0,0.0),Point_2(0.0,1.0),Point_2(-1.0,1.0)));

  Line_2 line;
  FT i = linear_least_squares_fitting_2(triangles.begin(),triangles.end(),line);

  std::cout<<"accuracy: "<<i<<std::endl;
  std::cout<<line<<std::endl;
  return 0;
}

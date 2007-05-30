// Example program for the linear_least_square_fitting function

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line_2;
typedef K::Point_2           Point_2;
typedef K::Circle_2          Circle_2;

int main()
{
  std::list<Circle_2> circles;
  circles.push_back(Circle_2(Point_2(0.0,0.0),9));
  circles.push_back(Circle_2(Point_2(0.0,10.0),49));
  circles.push_back(Circle_2(Point_2(10.0,0.0),49));

  Line_2 line;
  FT i = linear_least_squares_fitting_2(circles.begin(),circles.end(),line,CGAL::PCA_dimension_1_tag());

  std::cout<<"accuracy: "<<i<<std::endl;
  std::cout<<line<<std::endl;
  return 0;
}

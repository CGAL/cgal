// Example program for the linear_least_square_fitting function

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>
typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line_2;
typedef K::Point_2           Point_2;

int main(int argc, char** argv)
{
  std::list<Point_2> points;
  points.push_back(Point_2(1.0,0.0));
  points.push_back(Point_2(0.0,1.0));
  points.push_back(Point_2(1.0,1.0));
  points.push_back(Point_2(-1.0,1.0));
  points.push_back(Point_2(-1.0,0.0));
  //  points.push_back(Point_2(0.0,1.0));

  Line_2 line;
  linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::PCA_dimension_0_tag());
  
  std::cout<<line<<std::endl;
  return 0;
}

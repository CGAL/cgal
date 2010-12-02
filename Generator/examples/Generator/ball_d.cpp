#include <iostream>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>

typedef CGAL::Cartesian_d<double>                           Kd;
typedef Kd::Point_d                                         Point;
typedef CGAL::Creator_uniform_d<std::vector<double>::iterator, Point>Creator_d;

int main ()
{
  int nb_points = 10;
  int dim =5;
  double size = 100.0;
  std::cout << "Generating "<<nb_points<<" random points in a"
	    <<" ball in "<<dim<<"D of center 0 and radius "<<size<<std::endl;
  std::vector<Point> v;
  v.reserve (nb_points);
  CGAL::Random_points_in_ball_d<Point> gen (dim, 100.0);
  for (int i = 0; i < nb_points; ++i)  v.push_back (*gen++);
  for (int i = 0; i < nb_points; ++i)  std::cout<<"     "<<v[i]<<std::endl;
  return 0;
}

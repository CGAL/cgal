#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>
#include <vector>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Line_2 Line_2;

typedef CGAL::Random_points_in_square_2<Point_2,
                                        CGAL::Creator_uniform_2< double, Point_2 > > Point_generator;


int main()
{
  Line_2 line(-4.2885603045067812e-18, 1, 250.73609999999996);

  Point_2 point(35.306000000000004, 250.69800000000001);
  std::cout.precision(17);
  std::cout << line.projection(point) << std::endl;

  std::vector<Point_2> points;
  
  Point_generator pg(1000);
  for(int i = 0; i < 100000000; i++){
    Point_2 p = *pg++;
    points.push_back(p);
    Point_2 pp = line.projection(p);
    if(CGAL::squared_distance(pp,line) > 0.01){
      std::cout << " projection of " << p << "  is " << pp << " with squared distance to line " << CGAL::squared_distance(pp,line) << std::endl;
    }
  }

  CGAL::Timer t;
  t.start();
  double x=0;
   for(int i = 0; i < 100000000; i++){
     Point_2 pp = line.projection(points[i]);
     x += pp.x();
   }
   t.stop();
   std::cout << x << std::endl;
   std::cout << t.time() << " sec." << std::endl;
  return 0;
}

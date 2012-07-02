#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Delaunay_hyperbolic_triangulation_2.h>
#include <CGAL/Triangulation_hyperbolic_traits_2.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_hyperbolic_traits_2<K> Gt;

typedef Gt::Point_2 Point_2;
typedef Gt::Circle_2 Circle_2;
typedef Gt::FT FT;

typedef CGAL::Delaunay_hyperbolic_triangulation_2<Gt> Dt;

int main()
{
  CGAL::Timer timer;
  typedef CGAL::Creator_uniform_2<FT, Point_2> Creator;
  
  FT r = 100;
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_disc(r);
  
  int n = 10000;
  std::vector<Point_2> pts;
  std::vector<Point_2>::iterator ip;
  
  // Generating n random points
  for (int i=0 ; i < n ; i++) {
    pts.at(i) = *in_disc;
    in_disc++;
  }
  
  timer.start();
  
  Dt dt = Dt(Gt(r));
  
  for(ip = pts.begin(); ip != pts.end(); ++ip) {
    dt.insert(*ip);
  }
  
  timer.stop();
  
  std::cout << "Number of points: " << n << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  timer.reset();
                                                      
  return 0;
}

#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_triangulation_traits_2.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Hyperbolic_triangulation_traits_2<K> Gt;

typedef Gt::Point_2 Point_2;
typedef Gt::Circle_2 Circle_2;
typedef Gt::FT FT;

typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt> Dt;

int main()
{
  CGAL::Timer timer;
  typedef CGAL::Creator_uniform_2<FT, Point_2> Creator;
  
  FT r = 100;
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_disc(r);
  
  int n = 10000;
  std::cout << "Number of points: " << n << std::endl;

  std::vector<Point_2> pts(n);
  std::vector<Point_2>::iterator ip;
  
  // Generating n random points
  for (int i=0 ; i < n ; i++) {
    pts.at(i) = *in_disc;
    in_disc++;
  }

  std::cout << "check for hyperbolic faces during insertion" << std::endl;

  timer.start();
  
  Dt dt_during = Dt(Gt(r));
  
  for(ip = pts.begin(); ip != pts.end(); ++ip) {
    dt_during.insert(*ip);
  }
  
  timer.stop();
  
  assert(dt_during.is_valid());

  std::cout << "Number of vertices: " << dt_during.number_of_vertices() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  timer.reset();

  timer.start();
  
  std::cout << "check for hyperbolic faces only at the end" << std::endl;

  Dt dt_end = Dt(Gt(r));
  
  dt_end.insert(pts.begin(),pts.end());
  
  timer.stop();
  
  assert(dt_end.is_valid());

  std::cout << "Number of vertices: " << dt_end.number_of_vertices() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  timer.reset();
  
                                                      
  return 0;
}

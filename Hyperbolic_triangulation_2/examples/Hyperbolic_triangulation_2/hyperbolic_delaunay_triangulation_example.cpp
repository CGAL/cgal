#include <fstream>

// CGAL headers
#include <CGAL/IO/io.h>

#include <CGAL/Exact_circular_kernel_2.h>

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_circular_kernel_2             K;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<K> Gt;

typedef Gt::Point_2 Point_2;
typedef Gt::Circle_2 Circle_2;
typedef Gt::FT FT;

typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt> Dt;

int main()
{
  CGAL::Timer timer;
  typedef CGAL::Creator_uniform_2<FT, Point_2> Creator;
  
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_disc(FT(1));
  
  int n = 100;
  std::cout << "Number of points: " << n << std::endl;

  std::vector<Point_2> pts(n);
  std::vector<Point_2>::iterator ip;
  
  // Generating n random points
  for (int i=0 ; i < n ; i++) {
    pts.at(i) = *in_disc;
    in_disc++;
  }

  // std::vector<Point_2> pts;
  // std::vector<Point_2>::iterator ip;
  // Point_2 p;

  // std::ifstream ifs("input-file");
  // while(ifs >> p) {
  //   pts.push_back(p);
  // }
  // std::cout << "number of points " << std::distance(pts.begin(),pts.end()) << std::endl << std::endl;

  std::cout << "check for hyperbolic faces during insertion" << std::endl;

  timer.start();
  
  Dt dt_during;
  
  for(ip = pts.begin(); ip != pts.end(); ++ip) {
    dt_during.insert(*ip);
  }
  
  timer.stop();
  
  assert(dt_during.is_valid());

  std::cout << "Number of (finite) vertices: " << dt_during.number_of_vertices() << std::endl;
  std::cout << "number of (finite) Euclidean faces: " << dt_during.number_of_faces() << std::endl;
  std::cout << "number of hyperbolic faces: " << dt_during.number_of_hyperbolic_faces() << std::endl;
  std::cout << "number of hyperbolic edges: " << dt_during.number_of_hyperbolic_edges() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl << std::endl;
  
  timer.reset();

  timer.start();
  
  std::cout << "check for hyperbolic faces only at the end" << std::endl;

  Dt dt_end;
  
  dt_end.insert(pts.begin(),pts.end());
  
  timer.stop();
  
  assert(dt_end.is_valid());

  std::cout << "Number of (finite) vertices: " << dt_end.number_of_vertices() << std::endl;
  std::cout << "number of (finite) Euclidean faces: " << dt_end.number_of_faces() << std::endl;
  std::cout << "number of hyperbolic faces: " << dt_end.number_of_hyperbolic_faces() << std::endl;
  std::cout << "number of hyperbolic edges: " << dt_end.number_of_hyperbolic_edges() << std::endl;
  std::cout << "Time: " << timer.time() << std::endl;
  
  timer.reset();
  
                                                      
  return 0;
}

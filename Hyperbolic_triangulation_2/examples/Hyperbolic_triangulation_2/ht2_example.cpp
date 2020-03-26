#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <vector>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  Gt;
typedef Gt::Point_2                                         Point_2;
typedef Gt::Circle_2                                        Circle_2;
typedef Gt::FT                                              FT;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<Gt>       Dt;
typedef CGAL::Creator_uniform_2<FT, Point_2>                Creator;

int main(int argc, char** argv)
{
  int N;
  if(argc < 2) {
    std::cout << "usage: " << argv[0] << " [number_of_points]" << std::endl;
    std::cout << "Defaulting to 100k points..." << std::endl;
    N = 100000;
  } else {
    N = atoi(argv[1]);
  }

  CGAL::Timer timer;
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_disc;
  std::vector<Point_2> pts;
  std::vector<Point_2>::iterator ip;

  for(int i=0; i<N; ++i)
    pts.push_back(*(in_disc++));

  Dt dt_during;
  std::cout << "Insertion of points one by one (hyperbolic filtering at each step)"  << std::endl;
  std::cout << "===================================================================" << std::endl;
  timer.start();
  for(ip = pts.begin(); ip != pts.end(); ++ip)
    dt_during.insert(*ip);

  timer.stop();

  std::cout << "Number of vertices:         " << dt_during.number_of_vertices() << std::endl;
  std::cout << "Number of hyperbolic faces: " << dt_during.number_of_hyperbolic_faces() << std::endl;
  std::cout << "Number of hyperbolic edges: " << dt_during.number_of_hyperbolic_edges() << std::endl;
  std::cout << "Time:                       " << timer.time() << std::endl << std::endl;

  Dt dt_end;
  std::cout << "Insertion of point set (hyperbolic filtering only once at the end)"  << std::endl;
  std::cout << "===================================================================" << std::endl;
  timer.reset();
  timer.start();
  dt_end.insert(pts.begin(),pts.end());
  timer.stop();

  std::cout << "Number of vertices:         " << dt_end.number_of_vertices() << std::endl;
  std::cout << "Number of hyperbolic faces: " << dt_end.number_of_hyperbolic_faces() << std::endl;
  std::cout << "Number of hyperbolic edges: " << dt_end.number_of_hyperbolic_edges() << std::endl;
  std::cout << "Time:                       " << timer.time() << std::endl;

  return EXIT_SUCCESS;
}

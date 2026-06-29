// #define CGAL_SMALL_UNORDERED_MAP_STATS
// #define CGAL_PROFILE
//#define CGAL_USE_SSE2_FABS
//#define CGAL_USE_SSE2_MAX
//#define CGAL_MSVC_USE_STD_FABS  // use this one with precise
#define CGAL_NDEBUG 1
#define NDEBUG 1


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>
#include <locale>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

typedef K::Point_3                                           Point_3;
typedef CGAL::Timer                                          Timer;


int main(int argc, char* argv[])
{
  std::locale loc = std::locale()
      .combine<std::numpunct<char>>(std::locale("en_US.UTF8"));
  std::cout.imbue(loc);

  int M = 10; // outer loop counter
  int Q = 10000; // number of consecutive queries


  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/ocean_r.xyz");
  if(argc > 2) {
    auto M_ = std::atoi(argv[2]);
    if(M_ <= 0) {
      std::cerr << "Invalid number of iterations: " << M_ << ". Using default value of " << M << std::endl;
    } else {
      M = M_;
    }
  }
  if(argc > 3) {
    auto Q_ = std::atoi(argv[3]);
    if(Q <= 0) {
      std::cerr << "Invalid number of iterations: " << Q_ << ". Using default value of " << Q << std::endl;
    } else {
      Q = Q_;
    }
  }
  std::ifstream in(filename.c_str());
  std::vector<Point_3> points;
  Point_3 p, q;

  while(in >> p ){
    points.push_back(p);
  }

  std::cout << points.size() << " points read\n";

  std::cout << "We run " << M << " times: All four consecutive points combined with the next " << Q << " points as query" << std::endl;
  Timer timer;
  timer.start();
  std::size_t count = 0, inside = 0;
  std::size_t bound = (std::min)(points.size()/2, std::size_t(Q));
  for(int k = 0; k < M; k++)
    for(int i = 0; i < points.size()/2; ++i)
      for(int j = i+4; j < i+4+bound ; ++j ){
        ++count;
        if(side_of_oriented_sphere(points[i+(j%4)],points[i+(1+j%4)], points[i+(2+j%4)], points[i+(3+j%4)], points[j]) == CGAL::ON_NEGATIVE_SIDE) ++inside;
      }

  std ::cout << inside << " inside of " << count << " calls" << std::endl;

  timer.stop();
  std::cout << "Time elapsed: " << timer.time() << " sec" << std::endl;
  return 0;
}

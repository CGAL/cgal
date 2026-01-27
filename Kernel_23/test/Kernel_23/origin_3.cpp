//#define CGAL_PROFILE
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <vector>
#include <iostream>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Timer Timer;

int main(int argc, char* argv[] )
{
  std::ifstream ifs((argc>1)? argv[1]:CGAL::data_file_path("points_3/cube.xyz"));

  std::vector<Point_3> points;
  Point_3 p;

  while(ifs >> p){
    points.push_back(p);
  }

  std::cout << "Orientation_3" << std::endl;
  Timer t;

  const std::size_t N = points.size()-3;

  const K::Orientation_3 orientation = K().orientation_3_object();

  int positive = 0;

  t.start();
  {
    std::cout << "overload with 4 points" << std::endl;
    for(std::size_t k = 0; k < 100; ++k)
      for(std::size_t i = 0; i < N; ++i){
       Point_3 o(CGAL::ORIGIN);
      if(orientation(o, points[i], points[i+1], points[i+2]) == CGAL::POSITIVE){
        ++positive;
      }
    }
  }
  t.stop();

  std::cout << t.time() << " sec." << std::endl;

  t.reset();
  t.start();
  {
    std::cout << "overload with origin and 3 points" << std::endl;
    for (std::size_t k = 0; k < 100; ++k)
      for(std::size_t i = 0; i < N; ++i){
      if(orientation(CGAL::ORIGIN, points[i], points[i+1], points[i+2]) == CGAL::POSITIVE){
        --positive;
     }
    }
  }
  t.stop();


  if(positive != 0){
    std::cout << "Not the same results for Orientation_3"<< std::endl;
    assert(false);
  }
  std::cout << t.time() << " sec." << std::endl;


  std::cout << "Construct_orthogonal_vector_3" << std::endl;

  const K::Construct_orthogonal_vector_3 construct_orthogonal_vector = K().construct_orthogonal_vector_3_object();

  double sumx1 = 0, sumx2 = 0;

 t.start();
  {
    std::cout << "overload with 3 points" << std::endl;
    for(std::size_t k = 0; k < 100; ++k)
      for(std::size_t i = 0; i < N; ++i){
        Point_3 o(CGAL::ORIGIN);
        Vector_3 v = construct_orthogonal_vector(o, points[i], points[i+1]);
        sumx1 += CGAL::to_double(v.approx().x());
      }
  }
  t.stop();

  std::cout << t.time() << " sec." << std::endl;

  t.reset();
  t.start();
  {
    std::cout << "overload with origin and 2 points" << std::endl;
    for (std::size_t k = 0; k < 100; ++k)
      for(std::size_t i = 0; i < N; ++i){
        Vector_3 v = construct_orthogonal_vector(CGAL::ORIGIN, points[i], points[i+1]);
        sumx2 += CGAL::to_double(v.approx().x());
      }
  }

  t.stop();

  if(sumx1 != sumx2){
    std::cout << "Not the same results for Construct_orthogonal_vector" << std::endl;
    assert(false);
  }
  std::cout << t.time() << " sec." << std::endl;


  return 0;
}

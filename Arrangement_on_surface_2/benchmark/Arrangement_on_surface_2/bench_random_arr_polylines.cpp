#include <list>
#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Traits_2;
typedef Traits_2::Point_2                                 Point_2;
typedef Segment_traits_2::Curve_2                         Segment_2;
typedef Traits_2::Curve_2                                 Polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                     Arrangement_2;

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <number of points> [seed]"
              << std::endl;
    return -1;
  }
  unsigned int number_of_points(boost::lexical_cast<unsigned int>(argv[1]));
  std::list<Point_2> pts;
  unsigned int seed;
  if (argc == 3) {
    seed = boost::lexical_cast<unsigned int>(argv[2]);
    CGAL::Random rnd(seed);
    CGAL::Random_points_in_square_2<Point_2> g(10, rnd);
    for (unsigned int i = 1; i < number_of_points; ++i) pts.push_back(*g++);
  }
  else {
    CGAL::Random rnd;
    seed = rnd.get_seed();
    CGAL::Random_points_in_square_2<Point_2> g(10, rnd);
    for (unsigned int i = 1; i < number_of_points; ++i) pts.push_back(*g++);
  }
  std::cout << "Seed to be used: " << seed << std::endl;
  Polyline_2 poly(pts.begin(), pts.end());
  Arrangement_2 arr;
  boost::timer timer;
  insert(arr, poly);
  double secs = timer.elapsed();

  std::cout << "Arrangement computation took: " << secs << std::endl;
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Output for org-mode table
  std::cout << "| | "  << " | "
            << number_of_points << " | "
            << arr.number_of_vertices() << " | "
            << arr.number_of_edges() << " | "
            << arr.number_of_faces() << " |"
            << secs << " | " << std::endl;


  return 0;
}

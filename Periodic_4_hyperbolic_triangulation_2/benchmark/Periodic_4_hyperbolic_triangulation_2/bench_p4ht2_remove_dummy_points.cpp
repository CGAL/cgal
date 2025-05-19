#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/determinant.h>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <iostream>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>               Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>                Triangulation;
typedef Triangulation::Face_iterator                                                Face_iterator;
typedef Triangulation::Vertex_handle                                                Vertex_handle;
typedef Triangulation::Point                                                        Point;
typedef Traits::Side_of_original_octagon                                            Side_of_original_octagon;
typedef CGAL::Cartesian<double>::Point_2                                            Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                              Creator;

int main(int argc, char** argv)
{
  int iter;
  if(argc < 2)
  {
    std::cout << "usage: " << argv[0] << " [number_of_iterations]" << std::endl;
    std::cout << "defaulting to 10 iterations..." << std::endl;
    iter = 10;
  } else {
    iter = atoi(argv[1]);
  }

  Side_of_original_octagon pred;

  int N = 500;
  int min = 2*N;
  int max = -1;
  double mean = 0.0;
  for(int j=0; j<iter; ++j)
  {
    std::vector<Point_double> v;
    CGAL::Random_points_in_disc_2<Point_double, Creator> g(0.85);

    Triangulation tr;
    assert(tr.is_valid(true));

    int cnt = 0;
    int idx = 0;
    do
    {
      Point_double pt = *(++g);
      if(pred(pt) != CGAL::ON_UNBOUNDED_SIDE)
      {
        tr.insert(Point(pt.x(), pt.y()));
        cnt++;
      }
    }
    while(tr.number_of_dummy_points() > 0 && idx < N-1);

    if(tr.number_of_dummy_points() > 0)
    {
      std::cout << "!!! FAILED to remove all dummy points after the insertion of " << N << " random points!" << std::endl;
      continue;
    }

    assert(tr.is_valid());
    std::cout << cnt << std::endl;

    if(cnt > max)
      max = cnt;

    if(cnt < min)
      min = cnt;

    mean += cnt;
  }
  mean /= double(iter);

  std::cout << "Finished " << iter << " iterations!" << std::endl;
  std::cout << "Minimum number of points inserted: " << min << std::endl;
  std::cout << "Maximum number of points inserted: " << max << std::endl;
  std::cout << "Average number of points inserted: " << mean << std::endl;

  return EXIT_SUCCESS;
}

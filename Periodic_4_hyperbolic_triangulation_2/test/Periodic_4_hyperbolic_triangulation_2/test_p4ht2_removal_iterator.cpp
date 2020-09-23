#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>

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
typedef Triangulation::Vertex_iterator                                              Iter;
typedef Traits::Side_of_original_octagon                                            Side_of_original_octagon;

typedef CGAL::Cartesian<double>                                                     DKernel;
typedef DKernel::Point_2                                                            DPoint;
typedef CGAL::Creator_uniform_2<double, DPoint >                                    Creator;

int main(int argc, char** argv)
{
  int N;
  if(argc < 2)
  {
    std::cout << "usage: " << argv[0] << " [number of points]" << std::endl;
    std::cout << "generating 1000 points (default)!" << std::endl;
    N = 1000;
  }
  else
  {
    N = atoi(argv[1]);
  }

  int iters = 1;
  if(argc == 3)
    iters = atoi(argv[2]);


  for(int itr = 0; itr < iters; ++itr)
  {
    Triangulation tr;

    CGAL::Random_points_in_disc_2<DPoint, Creator> g(0.85);
    Side_of_original_octagon pred;

    std::vector<Vertex_handle> vh;

    int cnt = 0;
    do {
      DPoint pt = *g;
      ++g;
      if(pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
        vh.push_back(tr.insert(Point(pt.x(), pt.y())));
        cnt++;
      }
    }
    while(cnt < N);

    tr.try_to_remove_dummy_vertices();
    assert(tr.is_valid());

    tr.remove(vh.begin(), vh.end());

    std::cout << "Final count of vertices: " << tr.number_of_vertices() << std::endl;
    assert(tr.is_valid());
  }

  return EXIT_SUCCESS;
}

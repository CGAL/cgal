#include <CGAL/Cartesian.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/determinant.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>           Traits;
typedef Traits::FT                                                              NT;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef CGAL::Hyperbolic_octagon_translation_matrix<NT>                         Octagon_matrix;
typedef Triangulation::Point                                                    Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_original_octagon                                        Side_of_original_octagon;

typedef CGAL::Cartesian<double>::Point_2                                        Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                          Creator;

typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt             EKernel;
typedef CGAL::Delaunay_triangulation_2<EKernel>                                 Euclidean_triangulation;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<EKernel>               Ptraits;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<Ptraits>                      PEuclidean_triangulation;

typedef double                                                                  dNT;
typedef CGAL::Cartesian<dNT>                                                    dKernel;
typedef CGAL::Delaunay_triangulation_2<dKernel>                                 dTriangulation;

int main(int argc, char** argv)
{
  int N, iters;
  iters = 1;
  if(argc < 2)
  {
    std::cout << "usage: " << argv[0] << " [number_of_points_to_insert] [optional: number_of_iterations]" << std::endl;
    std::cout << "Defaulting to values: 10000000, 10..." << std::endl;
    N = 1000000;
    iters = 10;
  } else {
    N = atoi(argv[1]);
    if (argc < 3)
      iters = 1;
    else
      iters = atoi(argv[2]);
  }


  Side_of_original_octagon pred;

  std::cout << "---- for best results, make sure that you have compiled me in Release mode ----" << std::endl;

  double extime1 = 0.0;
  double extime2 = 0.0;
  double extime3 = 0.0;

  for(int exec = 1; exec <= iters; ++exec)
  {
    std::vector<Point> pts;
    std::vector<Point_double> dpts;
    CGAL::Random_points_in_disc_2<Point_double, Creator> g(0.85);

    int cnt = 0;
    std::cout << "================ iteration " << exec << " : generating points ================" << std::endl;
    do
    {
      Point_double pd = *(++g);
      Point pt = Point(pd.x(), pd.y());
      if(pred(pt) != CGAL::ON_UNBOUNDED_SIDE)
      {
        pts.push_back(pt);
        dpts.push_back(pd);
        ++cnt;
      }
    }
    while(cnt < N);

    if(cnt < N)
    {
      std::cout << "Failed to generate all the random points! Exiting..." << std::endl;
      return -1;
    }

    std::cout << "Inserting into hyperbolic periodic  CORE  triangulation...    "; std::cout.flush();
    Triangulation tr;
    CGAL::Timer t1;
    t1.start();
    tr.insert(pts.begin(), pts.end());
    t1.stop();
    extime1 += t1.time();
    std::cout << "DONE! (# of vertices = " << tr.number_of_vertices() << ", time = " << t1.time() << " secs)" << std::endl;

    std::cout << "inserting into Euclidean non-periodic  CORE  triangulation... "; std::cout.flush();
    Euclidean_triangulation etr;
    CGAL::Timer t2;
    t2.start();
    etr.insert(pts.begin(), pts.end());
    t2.stop();
    extime2 += t2.time();
    std::cout << "DONE! (# of vertices = " << etr.number_of_vertices() << ", time = " << t2.time() << " secs)" << std::endl;

    std::cout << "Inserting into Euclidean non-periodic DOUBLE triangulation... "; std::cout.flush();
    dTriangulation dtr;
    CGAL::Timer t3;
    t3.start();
    dtr.insert(dpts.begin(), dpts.end());
    t3.stop();
    extime3 += t3.time();
    std::cout << "DONE! (# of vertices = " << dtr.number_of_vertices() << ", time = " << t3.time() << " secs)" << std::endl;
  }

  double diters(iters);
  extime1 /= diters;
  extime2 /= diters;
  extime3 /= diters;

  std::cout << "Hyperbolic periodic      CORE  triangulation: average time = " << extime1 << std::endl;
  std::cout << "Euclidean  non-periodic  CORE  triangulation: average time = " << extime2 << std::endl;
  std::cout << "Euclidean  non-periodic DOUBLE triangulation: average time = " << extime3 << std::endl;

  return EXIT_SUCCESS;
}

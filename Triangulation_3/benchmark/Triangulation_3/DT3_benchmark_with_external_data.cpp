//#define CGAL_DO_NOT_USE_BOOST_MP

#define CGAL_LAZY_FILTERED_RATIONAL_KERNEL
#define CGAL_USE_FILTERED_RATIONAL_KERNEL

//#define CGAL_PROFILE
//#define CGAL_USE_SSE2_FABS
//#define CGAL_USE_SSE2_MAX
//#define CGAL_MSVC_USE_STD_FABS  // use this one with precise

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>

#include <chrono>
#include <iostream>
#include <fstream>

// typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Exact_predicates_exact_constructions_kernel       K;
typedef K::Point_3                                              Point_3;

typedef CGAL::Delaunay_triangulation_3<K>                       DT;

typedef CGAL::Triangulation_data_structure_3<
          CGAL::Triangulation_vertex_base_3<K>,
          CGAL::Delaunay_triangulation_cell_base_3<K>,
          CGAL::Parallel_tag>                                   PTds;
typedef CGAL::Delaunay_triangulation_3<K, PTds>                 PDT;

typedef CGAL::Timer                                             Timer;
typedef CGAL::Real_timer                                        Real_timer;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "usage: " << argv[0] << " data.xyz [parallel]" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::get_default_random() = CGAL::Random(0);
  std::cout << "seed:  " << CGAL::get_default_random().get_seed() << std::endl;

  std::cout << "FT: " << typeid(K::FT).name() << std::endl;
  int n;
  double x,y,z;
  std::vector<Point_3> points;
  Point_3 p;

  std::ifstream in(argv[1]);
  in >> n;
  points.reserve(n);
  while(n--)
  {
    in >> x >> y >> z;
    points.emplace_back(x,y,z);
  }

  int V, C;
  Timer timer;
  Real_timer rtimer;
  timer.start();
  rtimer.start();

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  bool use_parallel = false;
  if(argc > 2)
  {
    std::stringstream ss(argv[2]);
    if(!(ss >> std::boolalpha >> use_parallel))
    {
      std::cout << "second argument not a Boolean?" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if(use_parallel)
  {
    std::cout << "Parallel" << std::endl;
    CGAL::Bbox_3 bb  = CGAL::bounding_box(points.begin(), points.end()).bbox();

    PDT::Lock_data_structure locking_ds(bb, 50);

    PDT pdt(points.begin(), points.end(), &locking_ds);
    V = pdt.number_of_vertices();
    C = pdt.number_of_cells();
  }
  else
  {
    std::cout << "Sequential" << std::endl;
    DT dt(points.begin(), points.end());
    V = dt.number_of_vertices();
    C = dt.number_of_cells();
  }

  timer.stop();
  rtimer.stop();

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
          << "ms" << std::endl;

  std::cerr << "|V| = " <<  V << " |C| = " << C  << std::endl
            << timer.time() << " sec" << std::endl
            << rtimer.time() << " sec (r)" << std::endl;

  return EXIT_SUCCESS;
}

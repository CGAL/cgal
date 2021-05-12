#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/IO/Triangulation_off_ostream.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
#include <CGAL/Memory_sizer.h>

#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

//#define USE_DYNAMIC_KERNEL
#define OUTPUT_STATS_IN_CSV
//#define EXPORT_POINTS_TO_A_FILE

#ifdef OUTPUT_STATS_IN_CSV
static std::ofstream csv_file("stats.csv");
#endif

// Return the number of Bytes used
template<int D>
std::size_t compute_triangulation(std::size_t N)
{
#ifdef USE_DYNAMIC_KERNEL
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
#else
  typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > K;
#endif
  typedef CGAL::Delaunay_triangulation<K> DT;

  typedef typename DT::Vertex Vertex;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Full_cell Full_cell;
  typedef typename DT::Full_cell_handle Full_cell_handle;
  typedef typename DT::Facet Facet;
  typedef typename DT::Point Point;
  typedef typename DT::Geom_traits::RT RT;
  typedef typename DT::Finite_full_cell_const_iterator Finite_full_cell_const_iterator;

  typedef CGAL::Random_points_in_cube_d<Point> Random_points_iterator;
  CGAL::Timer cost;  // timer

  // Generate points
  std::vector<Point> points;
  CGAL::Random rng;
  Random_points_iterator rand_it(D, 2.0, rng);
  std::copy_n(rand_it, N, std::back_inserter(points));

#ifdef EXPORT_POINTS_TO_A_FILE
  std::ofstream os("points.txt");
  for (auto const& p : points)
  {
    CGAL::Triangulation_IO::output_point(os, K(), p);
    os << std::endl;
  }
#endif

  std::size_t mem_before = CGAL::Memory_sizer().virtual_size();
  cost.reset();
  cost.start();

  std::cout << "Delaunay triangulation of " << N <<
    " points in dim " << D << ":" << std::endl;

  DT dt(D);
  dt.insert(points.begin(), points.end());

  std::size_t mem = CGAL::Memory_sizer().virtual_size() - mem_before;
  double timing = cost.time();
  std::cout << "  Done in " << timing << " seconds." << std::endl;
  std::cout << "  Memory consumption: " << (mem >> 10) << " KB.\n";
  std::size_t nbfc= dt.number_of_finite_full_cells();
  std::size_t nbc= dt.number_of_full_cells();
  std::cout << "  " << dt.number_of_vertices() << " vertices, "
            << nbfc << " finite simplices and "
            << (nbc-nbfc) << " convex hull Facets."
            << std::endl;

#ifdef OUTPUT_STATS_IN_CSV
  csv_file
    << D << ";"
    << N << ";"
    << timing << ";"
    << mem << ";"
    << nbfc << "\n"
    << std::flush;
#endif


  return mem;
}

// Will compute triangulations of i*num_points_steps points,
// with i in [1, 2...], stopping after the last computation that takes
// more memory than mem_threshold_in_bytes
template<int D>
void go(
  std::size_t num_points_increment,
  std::size_t mem_threshold_in_MB = (3 << 10)) // 3 GB
{
  std::size_t mem = 0;
  for (std::size_t i = 1 ; mem < (mem_threshold_in_MB << 20) ; ++i)
  {
    mem = compute_triangulation<D>(i*num_points_increment);
  }
}

int main(int argc, char **argv)
{
  srand(static_cast<unsigned int>(time(NULL)));
  //int N = 100; if( argc > 1 ) N = atoi(argv[1]);
  go<2>(5000000);
  //go<3>(1000000);
  //go<4>(300000);
  //go<5>(50000);
  //go<6>(5000);
  //go<7>(1000);
  //go<8>(300);
  //go<9>(100);
  //go<10>(30);
  //go<11>(20);
  //go<12>(15);

  return 0;
}

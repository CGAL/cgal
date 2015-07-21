#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <string>


template<typename DTd, typename DT_23, typename RPI_23>
void test(int d, int N, std::string const& DTd_static_or_dyn)
{
  unsigned int seed = static_cast<unsigned int>(time(NULL));
  std::cerr << "Delaunay triangulation of " << N << " points in dim " << d;

  // Td
  {
    typedef typename DTd::Point Point;
    typedef CGAL::Random_points_in_cube_d<Point> Random_points_iterator;

    CGAL::default_random = CGAL::Random(seed);
    std::vector<Point> points;
    Random_points_iterator rand_it(d, 2.0);
    CGAL::cpp11::copy_n(rand_it, N, std::back_inserter(points));

    DTd dt(d);
    CGAL::Timer timer;
    timer.start();
    dt.insert(points.begin(), points.end());

    std::cerr << "  * Td: " << timer.time() << " s." << std::endl;
    std::size_t nbfc = dt.number_of_finite_full_cells();
    std::size_t nbc = dt.number_of_full_cells();
    std::cerr << "    " << dt.number_of_vertices() << " vertices, "
      << nbfc << " finite simplices and "
      << (nbc - nbfc) << " convex hull Facets."
      << std::endl;
  }

  // T2 or T3
  {
    typedef typename DT_23::Point Point_23;
    CGAL::default_random = CGAL::Random(seed); // Same seed
    std::vector<Point_23> points;
    RPI_23 rand_it(2.0);
    CGAL::cpp11::copy_n(rand_it, N, std::back_inserter(points));
    
    CGAL::Timer timer;
    timer.start();

    DT_23 dt;
    dt.insert(points.begin(), points.end());

    std::cerr << "  * T" << d << ": " << timer.time() << " s." << std::endl;
    /*std::size_t nbfc = dt.number_of_finite_full_cells();
    std::size_t nbc = dt.number_of_full_cells();
    std::cerr << "    " << dt.number_of_vertices() << " vertices, "
      << nbfc << " finite simplices and "
      << (nbc - nbfc) << " convex hull Facets."
      << std::endl;*/
  }

}

template< int D >
void go(const int N)
{
  CGAL_assertion(D == 2 || D == 3);

  // Td
  //typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > Kd;
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kd;
  typedef CGAL::Delaunay_triangulation<Kd> DTd;

  // T2 or T3
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K23;
  if (D == 2)
  {
    typedef CGAL::Delaunay_triangulation_2<K23> DT_2;
    typedef typename DT_2::Point Point;
    typedef CGAL::Random_points_in_square_2<Point> RPI_2;
    test<DTd, DT_2, RPI_2>(D, N, "static");
  }
  else if (D == 3)
  {
    typedef CGAL::Delaunay_triangulation_3<K23> DT_3;
    typedef typename DT_3::Point Point;
    typedef CGAL::Random_points_in_cube_3<Point> RPI_3;
    test<DTd, DT_3, RPI_3>(D, N, "static");
  }
}

int main(int argc, char **argv)
{
  srand(static_cast<unsigned int>(time(NULL)));
#ifdef _DEBUG
  int N = 100;
#else
  int N = 100000;
#endif
  if (argc > 1) N = atoi(argv[1]);
  go<2>(N);
  go<3>(N);

  return 0;
}

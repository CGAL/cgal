#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/IO/Triangulation_off_ostream.h>
#include <CGAL/algorithm.h>

#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

template<typename K, typename RTri>
void test(const int d, const string & type, const int N)
{
  typedef typename RTri::Vertex_handle Vertex_handle;
  typedef typename K::Weighted_point_d Weighted_point;
  typedef typename K::Point_d Bare_point;

  RTri rt(d);
  RTri rt_star_only(d);
  cerr << "\nBuilding Regular triangulation of (" << type << d
       << ") dimension with " << N << " points\n";
  assert(rt.empty());
  assert(rt_star_only.empty());

  srand(static_cast<unsigned int>(time(nullptr)));

  // Insert first point (0, 0...)
  vector<double> coords(d);
  for( int j = 0; j < d; ++j )
    coords[j] = 0;

  Weighted_point p = Weighted_point(
    Bare_point(d, coords.begin(), coords.end()),
    static_cast<double>(rand() % 10000)/100000);

  rt.insert(p);
  Vertex_handle first_vertex = rt_star_only.insert(p);

  // Insert the other points
  for( int i = 1 ; i < N ; ++i )
  {
    for( int j = 0; j < d; ++j )
      coords[j] = 10.*(rand() % RAND_MAX)/RAND_MAX - 5.;

    p = Weighted_point(
      Bare_point(d, coords.begin(), coords.end()),
      static_cast<double>(rand() % 10000)/1000000);

    rt.insert(p);
    rt_star_only.insert_if_in_star(p, first_vertex);
  }

  cerr << "\nChecking topology and geometry..."
       << (rt.is_valid(true) ? "OK.\n" : "Error.\n");

  cerr << "\nThe triangulation using 'insert' has current dimension " << rt.current_dimension()
       << " and " << rt.number_of_full_cells() << " full cells\n";

  cerr << "\nThe triangulation using 'insert_if_in_star' has current dimension " << rt.current_dimension()
       << " and " << rt_star_only.number_of_full_cells() << " full cells\n";

  // Export
  if (d <= 3)
  {
    std::ofstream off_stream_all("data/test_insert_all.off");
    CGAL::export_triangulation_to_off(off_stream_all, rt);
    std::ofstream off_stream_star_only("data/test_insert_if_in_star.off");
    CGAL::export_triangulation_to_off(off_stream_star_only, rt_star_only);
  }
}

template< int D >
void go(const int N)
{
  //typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> FK;
  typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > FK;
  typedef CGAL::Regular_triangulation<FK> Triangulation;
  //test<Triangulation>(D, "dynamic", N);
  test<FK, Triangulation>(D, "static", N);
}

int main()
{
  go<2>(100);
  return 0;
}

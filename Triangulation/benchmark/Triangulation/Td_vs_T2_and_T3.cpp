// To deactivate statics filters in the 2D/3D case
//#define CGAL_NO_STATIC_FILTERS

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Regular_triangulation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>

#include <vector>
#include <string>
#include "console_color.h"

template <typename DT_>
struct Stats_getter;

// T2 specialization
template <typename K>
struct Stats_getter<CGAL::Delaunay_triangulation_2<K> >
{
  typedef CGAL::Delaunay_triangulation_2<K> DT;

  Stats_getter(DT const& dt) : m_dt(dt) {}

  std::size_t number_of_vertices() { return m_dt.number_of_vertices(); }
  std::size_t number_of_finite_cells() { return m_dt.number_of_faces(); }

  DT m_dt;
};

// RT2 specialization
template <typename K>
struct Stats_getter<CGAL::Regular_triangulation_2<K> >
{
  typedef CGAL::Regular_triangulation_2<K> DT;

  Stats_getter(DT const& dt) : m_dt(dt) {}

  std::size_t number_of_vertices() { return m_dt.number_of_vertices(); }
  std::size_t number_of_finite_cells() { return m_dt.number_of_faces(); }

  DT m_dt;
};

// T3 specialization
template <typename K>
struct Stats_getter<CGAL::Delaunay_triangulation_3<K> >
{
  typedef CGAL::Delaunay_triangulation_3<K> DT;

  Stats_getter(DT const& dt) : m_dt(dt) {}

  std::size_t number_of_vertices() { return m_dt.number_of_vertices(); }
  std::size_t number_of_finite_cells() { return m_dt.number_of_finite_cells(); }

  DT m_dt;
};

// RT3 specialization
template <typename K>
struct Stats_getter<CGAL::Regular_triangulation_3<K> >
{
  typedef CGAL::Regular_triangulation_3<K> DT;

  Stats_getter(DT const& dt) : m_dt(dt) {}

  std::size_t number_of_vertices() { return m_dt.number_of_vertices(); }
  std::size_t number_of_finite_cells() { return m_dt.number_of_finite_cells(); }

  DT m_dt;
};


template<typename DT_d, typename DT_23,
         typename Pt_d_range, typename Pt_23_range>
void test(
  int d, int N, Pt_d_range const& points_d, Pt_23_range const& points_23,
  std::string const& DTd_static_or_dyn)
{
  // Td
  {
    DT_d dt(d);
    CGAL::Timer timer;
    timer.start();
    dt.insert(points_d.begin(), points_d.end());

    std::cerr << "  * Td: " << yellow << timer.time() << " s"
      << white << std::endl;
    std::cerr << "    " << dt.number_of_vertices() << " vertices, "
      << dt.number_of_finite_full_cells() << " finite cells."
      << std::endl;
  }

  // T2 or T3
  {
    CGAL::Timer timer;
    timer.start();

    DT_23 dt;
    dt.insert(points_23.begin(), points_23.end());

    std::cerr << "  * T" << d << ": " << yellow << timer.time() << " s"
      << white << std::endl;
    Stats_getter<DT_23> sg(dt);
    std::cerr << "    " << sg.number_of_vertices() << " vertices, "
      << sg.number_of_finite_cells() << " finite cells."
      << std::endl;
  }
}

template< int D, typename Dim_tag >
void go(const int N)
{
  CGAL_assertion(D == 2 || D == 3);

  // Generate points (in a common "array" format)
  std::vector<std::array<double, D> > coords;
  coords.reserve(N);
  for (int i = 0; i < N; ++i)
  {
    std::array<double, D> pt;
    for (int j = 0; j < D; ++j)
      pt[j] = CGAL::get_default_random().get_double(-1., 1.);
    coords.push_back(pt);
  }
  // Generate weights
  std::vector<double> weights;
  weights.reserve(N);
  for (int i = 0; i < N; ++i)
    weights.push_back(CGAL::get_default_random().get_double(-10., 10.));

  // DTd
  typedef CGAL::Epick_d<Dim_tag> Kd;
  typedef CGAL::Delaunay_triangulation<Kd> DT_d;
  typedef typename DT_d::Point Point_d;

  std::vector<Point_d> points_d;
  points_d.reserve(N);
  for (int i = 0; i < N; ++i)
    points_d.push_back(Point_d(D, coords[i].begin(), coords[i].end()));

  // RTd
  typedef CGAL::Regular_triangulation<Kd> RT_d;
  //typedef typename RT_d::Point Bare_point_d; // because of Regular_traits_adapter Point is actually a Weighted_point
  typedef typename RT_d::Weighted_point WPoint_d;

  std::vector<WPoint_d> wpoints_d;
  wpoints_d.reserve(N);
  for (int i = 0; i < N; ++i)
  {
    wpoints_d.push_back(WPoint_d(
      Point_d(D, coords[i].begin(), coords[i].end()),
      weights[i]));
  }

  // T2 or T3
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K23;
  if (D == 2)
  {
    // Delaunay
    typedef CGAL::Delaunay_triangulation_2<K23> DT_2;
    typedef typename DT_2::Point Point;

    std::vector<Point> points;
    points.reserve(N);
    for (int i = 0; i < N; ++i)
      points.push_back(Point(coords[i][0], coords[i][1]));

    std::cerr << std::endl << "DELAUNAY - dim " << D << " - "
      << N << " points." << std::endl;
    test<DT_d, DT_2>(D, N, points_d, points, "static");

    // Regular
    typedef CGAL::Regular_triangulation_2<K23> RT_2;
    typedef typename RT_2::Bare_point Bare_point;
    typedef typename RT_2::Point WPoint;

    std::vector<WPoint> wpoints;
    wpoints.reserve(N);
    for (int i = 0; i < N; ++i)
    {
      wpoints.push_back(WPoint(
        Bare_point(coords[i][0], coords[i][1]),
        weights[i]));
    }

    std::cerr << std::endl << "REGULAR - dim " << D << " - "
      << N << " points." << std::endl;
    test<RT_d, RT_2>(D, N, wpoints_d, wpoints, "static");
  }
  else if (D == 3)
  {
    typedef CGAL::Delaunay_triangulation_3<K23> DT_3;
    typedef typename DT_3::Point Point;

    std::vector<Point> points;
    points.reserve(N);
    for (int i = 0; i < N; ++i)
      points.push_back(Point(coords[i][0], coords[i][1], coords[i][2]));

    std::cerr << std::endl << "DELAUNAY - dim " << D << " - "
      << N << " points." << std::endl;
    test<DT_d, DT_3>(D, N, points_d, points, "static");

    // Regular
    typedef CGAL::Regular_triangulation_3<K23> RT_3;
    typedef typename RT_3::Bare_point Bare_point;
    typedef typename RT_3::Point WPoint;

    std::vector<WPoint> wpoints;
    wpoints.reserve(N);
    for (int i = 0; i < N; ++i)
    {
      wpoints.push_back(WPoint(
        Bare_point(coords[i][0], coords[i][1], coords[i][2]),
        weights[i]));
    }

    std::cerr << std::endl << "REGULAR - dim " << D << " - "
      << N << " points." << std::endl;
    test<RT_d, RT_3>(D, N, wpoints_d, wpoints, "static");
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

  std::cerr << "-----------------------------------------" << std::endl;
  std::cerr << "--             STATIC                  --" << std::endl;
  std::cerr << "-----------------------------------------" << std::endl;
  go<2, CGAL::Dimension_tag<2> >(N);
  go<3, CGAL::Dimension_tag<3> >(N);
  std::cerr << std::endl;

  std::cerr << "-----------------------------------------" << std::endl;
  std::cerr << "--             DYNAMIC                 --" << std::endl;
  std::cerr << "-----------------------------------------" << std::endl;
  go<2, CGAL::Dynamic_dimension_tag>(N);
  go<3, CGAL::Dynamic_dimension_tag>(N);
  std::cerr << std::endl;

  return 0;
}

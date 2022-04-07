// Synopsis: Testsuite of package Approximate_min_ellipsoid_d
//
// Revision: $Id$
// $Date$
//
// Author: Kaspar Fischer <fischerk@inf.ethz.ch> (ETH Zurich)

#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/algorithm.h>

#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_2.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_3.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>

template <typename T>
std::string tostr(const T& t)
{
  std::stringstream strm;
  strm << t;
  return strm.str();
}

template <typename T>
T fromstr(const std::string& s)
{
  std::stringstream strm(s);
  T t;
  strm >> t;
  return t;
}

typedef CGAL::MP_Float                                      ET;
typedef CGAL::Cartesian<double>                             K_2;
typedef CGAL::Approximate_min_ellipsoid_d_traits_2<K_2, ET> T_2;

typedef CGAL::Cartesian<double>                             K_3;
typedef CGAL::Approximate_min_ellipsoid_d_traits_3<K_3, ET> T_3;

typedef CGAL::Cartesian_d<double>                           K_d;
typedef CGAL::Approximate_min_ellipsoid_d_traits_d<K_d, ET> T_d;

void check(bool expr, const std::string& msg)
{
  if (!expr) {
    std::cerr << "\n" << msg << ": failed\n";
    std::exit(-1);
  }
}

struct TwoD {};
struct ThreeD {};
struct DD {};

template<typename Kernel,typename Point_list>
void add_random_points(int n, int d, int multiplicity, Point_list& list,
                       const DD&)
  // Adds n random d-dimensional points to the given list.
{
  typedef typename Point_list::value_type Point;
  CGAL::Random_points_in_cube_d<Point> rpg(d,100.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < multiplicity; ++j)
      list.push_back(*rpg);

    ++rpg;
  }
}

template<typename Kernel,typename Point_list>
void add_random_points(int n, int d, int multiplicity, Point_list& list,
                       const TwoD&)
  // Adds n random d-dimensional points to the given list.
{
  typedef typename Point_list::value_type Point;
  typedef CGAL::Point_d<K_d> Point_d;
  CGAL::Random_points_in_cube_d<Point_d> rpg(d,100.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < multiplicity; ++j) {
      Point_d p = *rpg;
      list.push_back(Point(*p.cartesian_begin(),*(p.cartesian_begin()+1)));
    }

    ++rpg;
  }
}

template<typename Kernel,typename Point_list>
void add_random_points(int n, int d, int multiplicity, Point_list& list,
                       const ThreeD&)
  // Adds n random d-dimensional points to the given list.
{
  typedef typename Point_list::value_type Point;
  typedef CGAL::Point_d<K_d> Point_d;
  CGAL::Random_points_in_cube_d<Point_d> rpg(d,100.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < multiplicity; ++j) {
      Point_d p = *rpg;
      list.push_back(Point(*p.cartesian_begin(),
                           *(p.cartesian_begin()+1),
                           *(p.cartesian_begin()+2)));
    }

    ++rpg;
  }
}

template<typename Kernel>
struct is_d_dimensional
{
  typedef DD value;
};

template<>
struct is_d_dimensional<T_2>
{
  typedef TwoD value;
};

template<>
struct is_d_dimensional<T_3>
{
  typedef ThreeD value;
};

int id = 0; // used to count files, see simple_test() below

template<typename Kernel, typename Traits>
void simple_test(int n, int d, int multiplicity, double eps)
  // Computes (1+eps)-approximations of MEL(P) for n random points in R^d.
{
  typedef typename Traits::Point                    Point;
  typedef std::vector<Point>                        Point_list;
  typedef CGAL::Approximate_min_ellipsoid_d<Traits> Mel;

  std::cerr << "n=" << n << ", d=" << d << ", mult=" << multiplicity
            << ", eps=" << eps;

  // generate points
  Point_list P;
  typedef typename is_d_dimensional<Traits>::value is_d;
  add_random_points<Kernel>(n, d, multiplicity, P, is_d());
  CGAL::cpp98::random_shuffle(P.begin(), P.end());

  // compute minellipsoid:
  Traits tco;
  Mel mel(eps, P.begin(), P.end(), tco);

  // check validity
  check(mel.is_valid(true), "validity");

  // find center:
  if (mel.is_full_dimensional()) {
    mel.center_cartesian_begin(); // (Note: forces center to be computed.)
    if (d == 2 || d == 3)
      mel.axes_lengths_begin();   // (Note: forces axes to be computed.)
  }

  // query
  bool is_fd   = mel.is_full_dimensional();
  check(!mel.is_empty() || !is_fd, "empty but full-dimensional");
  if (is_fd) {
    // call all accessors
    mel.achieved_epsilon();
    for (int i=0; i<d; ++i) {
      mel.defining_vector(i);
      for (int j=0; j<d; ++j)
        mel.defining_matrix(i,j);
    }
    mel.defining_scalar();
    std::cerr << ", achieved_eps=" << mel.achieved_epsilon();
    if (mel.achieved_epsilon() > eps)
      std::cerr << " (desired eps not achieved)";

    // write EPS
    if (d == 2) {
      mel.write_eps(tostr(id)+".eps");
      std::cerr << " (file " << id << ")";
      ++id;
    }
  }
  std::cerr << std::endl;
}

void test(int n, int multiplicity)
{
  // 2d
  simple_test<K_d,T_d>(n, 2, multiplicity, 0.5);
  simple_test<K_d,T_d>(n, 2, multiplicity, 0.1);
  simple_test<K_d,T_d>(n, 2, multiplicity, 0.01);
  simple_test<K_d,T_d>(n, 2, multiplicity, 0.001);

  simple_test<K_2,T_2>(n, 2, multiplicity, 0.5);
  simple_test<K_2,T_2>(n, 2, multiplicity, 0.1);
  simple_test<K_2,T_2>(n, 2, multiplicity, 0.01);
  simple_test<K_2,T_2>(n, 2, multiplicity, 0.001);

  // 3d
  simple_test<K_d,T_d>(n, 3, multiplicity, 0.5);
  simple_test<K_d,T_d>(n, 3, multiplicity, 0.1);
  simple_test<K_d,T_d>(n, 3, multiplicity, 0.01);
  simple_test<K_d,T_d>(n, 3, multiplicity, 0.001);

  simple_test<K_3,T_3>(n, 3, multiplicity, 0.5);
  simple_test<K_3,T_3>(n, 3, multiplicity, 0.1);
  simple_test<K_3,T_3>(n, 3, multiplicity, 0.01);
  simple_test<K_3,T_3>(n, 3, multiplicity, 0.001);

  // 5d
  simple_test<K_d,T_d>(n, 5, multiplicity, 0.5);
  simple_test<K_d,T_d>(n, 5, multiplicity, 0.1);
  simple_test<K_d,T_d>(n, 5, multiplicity, 0.01);
}

int main()
{
  // multiplicity 1
  test(1,   1);
  test(2,   1);
  test(5, 1);
  test(100, 1);

  // multiplicity 2
  test(1,   2);
  test(2,   2);
  test(100, 2);

  // high multiplicity
  test(2,   2);
  test(10, 10);
  test(100,10);
}

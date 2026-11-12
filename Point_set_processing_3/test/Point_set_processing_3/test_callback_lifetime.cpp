// Regression test: a `callback` passed as a `std::function<bool(double)>`
// must stay alive during the whole run of the algorithm.
//
// The algorithms used to bind the parameter, returned by value by
// `get_parameter()`, to a `const std::function<bool(double)>&`, which
// dangled. Depending on the platform, the callback was then silently
// skipped (MSVC), called on freed memory (GCC), or flagged by ASan.
//
// The callback owns a heap-allocated payload (too large for the small
// buffer optimization of std::function) and checks it at each call.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>

#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/cluster_point_set.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/jet_smooth_point_set.h>
#endif

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <vector>

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3   = Kernel::Point_3;
using Vector_3  = Kernel::Vector_3;
using Point_set = CGAL::Point_set_3<Point_3>;
using Callback  = std::function<bool(double)>;
using Tag       = CGAL::Sequential_tag;

bool all_ok = true;

// Calls `algorithm(points, callback)` on a noisy sphere and checks that the
// callback was called, with a valid state and a progress in [0,1].
template <typename Algorithm>
void test(const char* name, Algorithm algorithm)
{
  CGAL::Random rnd(42);
  Point_set points(true);
  for (int i = 0; i < 1000; ++i)
  {
    Vector_3 v(rnd.get_double(-1., 1.), rnd.get_double(-1., 1.), rnd.get_double(-1., 1.));
    v = v / std::sqrt(v.squared_length());
    points.insert(CGAL::ORIGIN + v * (1. + rnd.get_double(-0.01, 0.01)), v);
  }

  std::size_t nb_calls = 0;
  bool valid = true;
  const std::vector<double> payload(256, 42.);
  const Callback callback = [payload, &nb_calls, &valid](double advancement)
  {
    valid = valid && payload.size() == 256 && payload.back() == 42.
                  && advancement >= 0. && advancement <= 1.;
    ++nb_calls;
    return true;
  };

  algorithm(points, callback);

  const bool ok = nb_calls > 0 && valid;
  std::cout << (ok ? "[OK]   " : "[FAIL] ") << name << " (" << nb_calls << " calls"
            << (valid ? "" : ", invalid callback state or progress") << ")" << std::endl;
  all_ok = all_ok && ok;
}

int main()
{
  const unsigned int k = 12;

  test("compute_average_spacing", [&](Point_set& ps, const Callback& cb)
       { CGAL::compute_average_spacing<Tag>(ps, k, ps.parameters().callback(cb)); });
  test("remove_outliers", [&](Point_set& ps, const Callback& cb)
       { CGAL::remove_outliers<Tag>(ps, k, ps.parameters().callback(cb)); });
  test("pca_estimate_normals", [&](Point_set& ps, const Callback& cb)
       { CGAL::pca_estimate_normals<Tag>(ps, k, ps.parameters().callback(cb)); });
#ifdef CGAL_EIGEN3_ENABLED
  test("jet_estimate_normals", [&](Point_set& ps, const Callback& cb)
       { CGAL::jet_estimate_normals<Tag>(ps, k, ps.parameters().callback(cb)); });
  test("jet_smooth_point_set", [&](Point_set& ps, const Callback& cb)
       { CGAL::jet_smooth_point_set<Tag>(ps, k, ps.parameters().callback(cb)); });
#endif
  test("bilateral_smooth_point_set", [&](Point_set& ps, const Callback& cb)
       { CGAL::bilateral_smooth_point_set<Tag>(ps, k, ps.parameters().callback(cb)); });
  test("wlop_simplify_and_regularize_point_set", [&](Point_set& ps, const Callback& cb)
       {
         std::vector<Point_3> output;
         CGAL::wlop_simplify_and_regularize_point_set<Tag>(
           ps, std::back_inserter(output), ps.parameters().number_of_iterations(2).callback(cb));
       });
  test("hierarchy_simplify_point_set", [&](Point_set& ps, const Callback& cb)
       { CGAL::hierarchy_simplify_point_set(ps, ps.parameters().callback(cb)); });
  test("cluster_point_set", [&](Point_set& ps, const Callback& cb)
       {
         std::vector<std::pair<std::size_t, std::size_t> > adjacencies;
         CGAL::cluster_point_set(ps, ps.add_property_map<int>("cluster", -1).first,
                                 ps.parameters().neighbor_radius(0.1)
                                                .adjacencies(std::back_inserter(adjacencies))
                                                .callback(cb));
       });

  if (!all_ok)
  {
    std::cerr << "Some tests FAILED" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

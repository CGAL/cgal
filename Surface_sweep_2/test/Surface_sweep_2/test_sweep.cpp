// examples/Pm_with_intersections/example4
// ---------------------------------------
#include <algorithm>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <CGAL/config.h>

#define CGAL_SEGMENT_TRAITS          1
#define CGAL_POLYLINE_TRAITS         11
#define CGAL_CONIC_TRAITS            21
#define CGAL_POLYCONIC_TRAITS        22

// Picking a default Traits class (this, with the
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_TRAITS
#endif

// Making sure test doesn't fail if CORE is not installed
#if ! defined(CGAL_USE_CORE) && \
  ((CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS) ||       \
   (CGAL_ARR_TEST_TRAITS == CGAL_POLYCONIC_TRAITS))

int main() {
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}

#elif ! defined(CGAL_USE_GMP) && \
  ((CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS) || \
   (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS))

int main() {
  std::cout << "A try to run test with GMP number type but GMP is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}

#else

// Choose traits

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYCONIC_TRAITS
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#else
#error No traits defined for test
#endif

#include <list>
#include <string>

#include <CGAL/Surface_sweep_2_algorithms.h>

#include "Compare_curves.h"

#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS
typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>         Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Rat_kernel::Segment_2                           Rat_segment_2;
typedef Rat_kernel::Circle_2                            Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Traits_2;
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYCURVE_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Rat_kernel::Segment_2                           Rat_segment_2;
typedef Rat_kernel::Circle_2                            Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Conic_traits;
typedef CGAL::Arr_polycurve_traits_2<Conic_traits>      Traits;

#endif

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

typedef std::list<Point_2>                              Points;
typedef std::list<Curve_2>                              Curves;
typedef std::list<X_monotone_curve_2>                   X_monotone_curves;

bool read_curves(std::ifstream& inp, Curves& curves, const Traits& traits);
bool read_xcurves(std::ifstream& inp, X_monotone_curves& xcurves,
                            const Traits& traits);
bool read_points(std::ifstream& inp, Points& points, const Traits& traits);
bool curves_identical(X_monotone_curves& list1, X_monotone_curves& list2);
bool points_identical(Points& list1, Points& list2);

// istream modifier skips chars until end of line.
std::istream& skip_until_eol(std::istream& in) {
  if (in.eof()) return in;
  char c;
  while (in.get(c) && (c != '\n'));
  return in;
}

// istream modifier that checks for OFF comments and removes them.
std::istream& skip_comment(std::istream& in) {
  char c;
  while ((in >> c) && (c == '#')) in >> skip_until_eol;
  in.putback(c);
  return in;
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Specify a file name " << std::endl;
    return -1;
  }

  std::ifstream inp(argv[1]);
  if (! inp.is_open()) {
    std::cerr << "Error: Cannot open file " << argv[1] << "!" << std::endl;
    return -1;
  }

  Traits tr;

  Curves curves;
  if (! read_curves(inp, curves, tr)) return -1;

  // Test subcurves w/o overlapping
  X_monotone_curves curves_no_overlap_out;
  CGAL::compute_subcurves(curves.begin(), curves.end(),
                          std::back_inserter(curves_no_overlap_out),
                          false, tr);


  X_monotone_curves curves_no_overlap;
  if (! read_xcurves(inp, curves_no_overlap, tr)) return -1;

  if (! compare_lists(curves_no_overlap_out, curves_no_overlap, tr)) {
    std::cerr << "Curves w/o overlapping do not match!\n";
    for (const auto& xcv : curves_no_overlap_out) std::cerr << xcv << std::endl;
    return -1;
  }

  // Test intersection points (with endpoints)
  Points points_with_ends_out;
  CGAL::compute_intersection_points(curves.begin(), curves.end(),
                                    std::back_inserter(points_with_ends_out),
                                    true, tr);

  Points points_with_ends;
  if (! read_points(inp, points_with_ends, tr)) return -1;

  if (! compare_lists(points_with_ends_out, points_with_ends, tr)) {
    std::cerr << "Endpoints do not match!\n";
    for (const auto& p : points_with_ends_out) std::cerr << p << std::endl;
    return -1;
  }

  // Test intersection points w/o end points
  Points points_without_ends_out;
  CGAL::compute_intersection_points(curves.begin(), curves.end(),
                                    std::back_inserter(points_without_ends_out),
                                    false, tr);

  Points points_without_ends;
  if (! read_points(inp, points_without_ends, tr)) return -1;

  if (! compare_lists(points_without_ends_out, points_without_ends, tr)) {
    std::cerr << "Intersection points do not match!\n";
    for (const auto& p : points_without_ends_out) std::cerr << p << std::endl;
    return -1;
  }

  // Test subcurves w/ overlapping
  X_monotone_curves curves_with_overlap_out;
  CGAL::compute_subcurves(curves.begin(), curves.end(),
                          std::back_inserter(curves_with_overlap_out),
                          true, tr);

  X_monotone_curves curves_with_overlap;
  if (! read_xcurves(inp, curves_with_overlap, tr)) return -1;

  if (! compare_lists(curves_with_overlap_out, curves_with_overlap, tr)) {
    std::cerr << "Curves w/ overlapping do not match!\n";
    for (const auto& xcv : curves_with_overlap_out)
      std::cerr << xcv << std::endl;
    return -1;
  }

  // Test the do_curves_intersecting method
  bool do_intersect_out =
    CGAL::do_curves_intersect(curves.begin(), curves.end());

  bool do_intersect = false;
  if ((points_without_ends.size() != 0) ||
      (curves_no_overlap_out.size() != curves_with_overlap_out.size()))
    do_intersect = true;

  if (do_intersect_out != do_intersect) {
    std::cerr << "Error: do_intersect()\n";
    return -1;
  }

  std::cout << "Passed\n";
  return 0;
}

#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS

bool read_curves(std::ifstream& inp, Curves& curves, const Traits&) {
  int count;
  inp >> skip_comment >> count;
  std::cout << "read_curves " << count << "\n";

  for (int i = 0; i < count; ++i) {
    NT x0, y0, x1, y1;
    inp >> skip_comment >> x0 >> y0 >> x1 >> y1;
    Point_2 p1(x0, y0);
    Point_2 p2(x1, y1);
    Curve_2 curve(p1, p2);
    curves.push_back(curve);
    std::cout << curve << "\n";
  }
  return true;
}

bool read_xcurves(std::ifstream& inp, X_monotone_curves& xcurves,
                  const Traits& traits)
{ return read_curves(inp, xcurves, traits); }

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS

template <typename Curves_, typename Ctr>
bool read_curves_(std::ifstream& inp, Curves_& curves, const Traits& traits,
                  const Ctr& ctr) {
  int count;
  inp >> skip_comment >> count;
  // std::cout << "read_curves " << count << "\n";
  for (int i = 0; i < count; ++i) {
    Points points;
    auto rc = read_points(inp, points, traits);
    if (! rc) return false;
    auto cv = ctr(points.begin(), points.end());
    // std::cout << cv << "\n";
    curves.push_back(cv);
  }
  return true;
}

/*! Read curves.
 */
bool read_curves(std::ifstream& inp, Curves& curves, const Traits& traits) {
  auto ctr_cv = traits.construct_curve_2_object();
  return read_curves_(inp, curves, traits, ctr_cv);
}

/*! Read x-monotone curves.
 */
bool read_xcurves(std::ifstream& inp, X_monotone_curves& xcurves,
                  const Traits& traits) {
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();
  return read_curves_(inp, xcurves, traits, ctr_xcv);
}

#else
#error No traits defined for test
#endif

bool read_points(std::ifstream& inp, Points& points, const Traits&) {
  int count;
  inp >> skip_comment >> count;

  // std::cout << "read_points " << count << "\n";
  for (int i = 0; i < count; i++) {
    NT x, y;
    inp >> skip_comment >> x >> y;
    Point_2 p(x, y);
    // std::cout << p << "\n";
    points.push_back(p);
  }
  return true;
}

#endif

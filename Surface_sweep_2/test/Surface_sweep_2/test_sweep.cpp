// examples/Pm_with_intersections/example4
// ---------------------------------------
#include <algorithm>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <CGAL/config.h>
#include <CGAL/Arrangement_2.h>
// #include <CGAL/draw_arrangement_2.h>

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
                                                        Traits;
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYCONIC_TRAITS
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

#if CGAL_ARR_TEST_TRAITS != CGAL_CONIC_TRAITS

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

#if CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS

bool read_curves(std::ifstream& inp, Curves& curves, const Traits&) {
  int count;
  inp >> skip_comment >> count;
  // std::cout << "read_curves " << count << "\n";

  for (int i = 0; i < count; ++i) {
    NT x0, y0, x1, y1;
    inp >> skip_comment >> x0 >> y0 >> x1 >> y1;
    Point_2 p1(x0, y0);
    Point_2 p2(x1, y1);
    Curve_2 curve(p1, p2);
    curves.push_back(curve);
    // std::cout << curve << "\n";
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

#elif CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS

void read_curve(std::ifstream& is, Curve_2& cv, const Traits& traits) {
  // Read a line from the input file.
  char one_line[128];
  auto ctr_curve_2 = traits.construct_curve_2_object();

  is >> skip_comment;
  is.getline(one_line, 128);
  std::string stringvalues(one_line);
  std::istringstream str_line(stringvalues, std::istringstream::in);

  // Get the arc type.
  // Supported types are: 'f' - Full ellipse (or circle).
  //                      'e' - Elliptic arc (or circular arc).
  //                      's' - Line segment.
  char type;
  bool is_circle = false;               // Is this a circle.
  Rat_circle_2 circle;
  Rational r, s, t, u, v, w;            // The conic coefficients.

  str_line >> type;

  // An ellipse (full ellipse or a partial ellipse):
  if (type == 'f' || type == 'F' || type == 'e' || type == 'E') {
    // Read the ellipse (using the format "a b x0 y0"):
    //
    //     x - x0   2      y - y0   2
    //  ( -------- )  + ( -------- )  = 1
    //       a               b
    //
    int a, b, x0, y0;

    str_line >> a >> b >> x0 >> y0;

    Rational a_sq = Rational(a*a);
    Rational b_sq = Rational(b*b);

    if (a == b) {
      is_circle = true;
      circle =
        Rat_circle_2(Rat_point_2(Rational(x0), Rational(y0)), Rational(a*b));
    }
    else {
      r = b_sq;
      s = a_sq;
      t = 0;
      u = Rational(-2*x0*b_sq);
      v = Rational(-2*y0*a_sq);
      w = Rational(x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq);
    }

    if (type == 'f' || type == 'F') {
      // Create a full ellipse (or circle).
      if (is_circle) cv = ctr_curve_2(circle);
      else cv = ctr_curve_2(r, s, t, u, v, w);
    }
    else {
      // Read the endpointd of the arc.
      int x1, y1, x2, y2;

      str_line >> x1 >> y1 >> x2 >> y2;

      Point_2 source = Point_2(Algebraic(x1), Algebraic(y1));
      Point_2 target = Point_2(Algebraic(x2), Algebraic(y2));

        // Create the arc. Note that it is always clockwise oriented.
      if (is_circle) cv = ctr_curve_2(circle, CGAL::CLOCKWISE, source, target);
      else cv = ctr_curve_2(r, s, t, u, v, w, CGAL::CLOCKWISE, source, target);
    }
  }
  else if (type == 's' || type == 'S') {
    // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
    int x1, y1, x2, y2;

    str_line >> x1 >> y1 >> x2 >> y2;

    // Create the segment.
    Rat_point_2 source = Rat_point_2 (Rational(x1), Rational(y1));
    Rat_point_2 target = Rat_point_2 (Rational(x2), Rational(y2));

    cv = ctr_curve_2(Rat_segment_2 (source, target));
  }

  // std::cout << cv << std::endl;
}

/*! Read curves.
 */
bool read_curves(std::ifstream& inp, Curves& curves, const Traits& traits) {
  // auto ctr_cv = traits.construct_curve_2_object();
  int count;
  inp >> skip_comment >> count;
  Curve_2 cv;
  char dummy[256];
  inp.getline(dummy, sizeof(dummy));
  for (int i = 0; i < count; ++i) {
    read_curve(inp, cv, traits);
    curves.push_back(cv);
  }
  return true;
}

#else
#error No traits defined for test
#endif

#if CGAL_ARR_TEST_TRAITS != CGAL_CONIC_TRAITS

// Test subcurves w/o overlapping
bool test_curves_no_overlap(std::ifstream& inp, Curves& /* curves */,
                            const X_monotone_curves& curves_no_overlap_out,
                            const Traits& tr) {
  X_monotone_curves curves_no_overlap;
  if (! read_xcurves(inp, curves_no_overlap, tr)) return false;

  if (! compare_lists(curves_no_overlap_out, curves_no_overlap, tr)) {
    std::cerr << "Error: Curves w/o overlapping do not match!\n";
    for (const auto& xcv : curves_no_overlap_out) std::cerr << xcv << std::endl;
    return false;
  }

  return true;
}

// Test subcurves w/ overlapping
bool test_curves_with_overlap(std::ifstream& inp, Curves& /* curves */,
                              const X_monotone_curves& curves_with_overlap_out,
                              const Traits& tr) {
  X_monotone_curves curves_with_overlap;
  if (! read_xcurves(inp, curves_with_overlap, tr)) return false;

  if (! compare_lists(curves_with_overlap_out, curves_with_overlap, tr)) {
    std::cerr << "Error: Curves w/ overlapping do not match!\n";
    for (const auto& xcv : curves_with_overlap_out)
      std::cerr << xcv << std::endl;
    return false;
  }

  return true;
}

// Test intersection points (with endpoints)
bool test_points_with_ends(std::ifstream& inp, Curves& /* curves */,
                           const Points& points_with_ends_out,
                           const Traits& tr) {
  Points points_with_ends;
  if (! read_points(inp, points_with_ends, tr)) return false;

  if (! compare_lists(points_with_ends_out, points_with_ends, tr)) {
    std::cerr << "Error: Endpoints do not match!\n";
    for (const auto& p : points_with_ends_out) std::cerr << p << std::endl;
    return false;
  }

  return true;
}

// Test intersection points w/o end points
bool test_points_no_ends(std::ifstream& inp, Curves& /* curves */,
                         const Points& points_no_ends_out,
                         const Traits& tr) {
  Points points_no_ends;
  if (! read_points(inp, points_no_ends, tr)) return false;

  if (! compare_lists(points_no_ends_out, points_no_ends, tr)) {
    std::cerr << "Error: Intersection points do not match!\n";
    for (const auto& p : points_no_ends_out) std::cerr << p << std::endl;
    return false;
  }

  return true;
}

#else

/*! Test the surface sweep with conic traits.
 */
bool test_conic(std::ifstream& inp, Curves& curves,
                const X_monotone_curves& curves_no_overlap_out,
                const Points& points_with_ends_out,
                const Points& points_no_ends_out,
                const Traits& traits) {
  auto ctr_bbox_2 = traits.construct_bbox_2_object();
  CGAL::Bbox_2 bbox = ctr_bbox_2(curves.front());
  for (auto it = std::next(curves.begin()); it != curves.end(); ++it)
    bbox = bbox + ctr_bbox_2(*it);

  // generate the string for the output
  std::stringstream out1;
  for (const auto& xcv : curves_no_overlap_out) out1 << xcv << "\n";

  // read the output from the file
  std::stringstream out2;
  char buf[1024];
  int count = 0;

  inp >> count;
  inp.getline(buf, 1024); // to get rid of the new line
  for (int i = 0; i < count; ++i) {
    inp.getline(buf, 1024);
    out2 << buf << "\n";
  }

  // std::cout << "Result: \n" << curves_no_overlap_out.size() << "\n";
  // for (const auto& xcv : curves_no_overlap_out)
  //   std::cout << xcv << "\n";

  std::string calculated = out1.str();
  std::string infile = out2.str();

  if (infile != calculated) {
    std::cerr << "Error\n";
    std::cerr << "\ncalculated:\n";
    std::cerr << calculated << std::endl;
    std::cerr << "\nin file:\n";
    std::cerr << infile << std::endl;
    std::cerr << "--"  << std::endl;
    return false;
  }

  std::size_t points_with_ends_size, points_no_ends_size;
  inp >> skip_comment >> points_with_ends_size >> points_no_ends_size;

  auto points_with_ends_out_size = points_with_ends_out.size();
  if (points_with_ends_size != points_with_ends_out_size ) {
    std::cerr << "Error: Number of endpoints do not match ("
              << points_with_ends_out_size << ", "
              << points_with_ends_size << ")\n";
    return false;
  }

  auto points_no_ends_out_size = points_no_ends_out.size();
  if (points_no_ends_size != points_no_ends_out_size) {
    std::cerr << "Error: Number of intersection points do not match ("
              << points_no_ends_out_size << ", "
              << points_no_ends_size << ")\n";
    return false;
  }

  return true;
}

#endif

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

  // {
  //   using Arrangement = CGAL::Arrangement_2<Traits>;
  //   Arrangement arr(&tr);
  //   CGAL::insert(arr, curves.begin(), curves.end());
  //   CGAL::draw(arr, "conics", true);
  // }

  X_monotone_curves curves_no_overlap_out;
  CGAL::compute_subcurves(curves.begin(), curves.end(),
                          std::back_inserter(curves_no_overlap_out),
                          false, tr);

  X_monotone_curves curves_with_overlap_out;
  CGAL::compute_subcurves(curves.begin(), curves.end(),
                          std::back_inserter(curves_with_overlap_out),
                          true, tr);

  Points points_with_ends_out;
  CGAL::compute_intersection_points(curves.begin(), curves.end(),
                                    std::back_inserter(points_with_ends_out),
                                    true, tr);

  Points points_no_ends_out;
  CGAL::compute_intersection_points(curves.begin(), curves.end(),
                                    std::back_inserter(points_no_ends_out),
                                    false, tr);

#if CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS
  if (! test_conic(inp, curves, curves_no_overlap_out,
                   points_with_ends_out, points_no_ends_out, tr))
    return -1;
#else
  if (! test_curves_no_overlap(inp, curves, curves_no_overlap_out, tr))
    return -1;
  if (! test_points_with_ends(inp, curves, points_with_ends_out, tr)) return -1;
  if (! test_points_no_ends(inp, curves, points_no_ends_out, tr)) return -1;
  if (! test_curves_with_overlap(inp, curves, curves_with_overlap_out, tr))
    return -1;
#endif

  // Test the do_curves_intersecting method
  bool do_intersect_out =
    CGAL::do_curves_intersect(curves.begin(), curves.end());
  bool do_intersect = ! points_no_ends_out.empty() ||
    (curves_no_overlap_out.size() != curves_with_overlap_out.size());
  if (do_intersect_out != do_intersect) {
    std::cerr << "Error: do_intersect()\n";
    return -1;
  }

  std::cout << "Passed\n";
  return 0;
}

#endif

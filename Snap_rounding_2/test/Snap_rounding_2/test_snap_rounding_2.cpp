#include <algorithm>
#include <cstdlib>

#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Hot_pixel_snap_rounding_traits_2.h>
#include <CGAL/Double_grid_snap_rounding_traits_2.h>
#include <CGAL/Float_grid_snap_rounding_traits_2.h>
#include <CGAL/Integer_grid_snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <CGAL/Polygon_2.h>

#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/Snap_rounding_traits_2.h>
#endif // CGAL_NO_DEPRECATED_CODE

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Segment_2 = Kernel::Segment_2;
using Point_2 = Kernel::Point_2;
using Polyline_2 = std::vector<Point_2>;
using Vector_2 = Kernel::Vector_2;
using FT = Kernel::FT;
using Float_grid_traits = CGAL::Float_grid_snap_rounding_traits_2<Kernel>;
using Hot_pixel_traits = CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>;

namespace params = CGAL::parameters;

void test_api(){
  std::vector< Segment_2 > segs;
  std::vector< boost::container::small_vector<Point_2, 2> > out_segs;
  std::vector< Polyline_2 > out_poly;

  FT e(std::pow(2, -60));
  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
  segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
  segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));

  // Default
  snap_rounding_2(segs, out_segs, params::output_unique_segments(true));
  snap_rounding_2(segs, out_poly);
  assert(out_segs.size() == 11);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();

  // VSSR Traits
  snap_rounding_2(segs, out_segs, params::output_unique_segments(true).geom_traits(Float_grid_traits()));
  snap_rounding_2(segs, out_poly, params::geom_traits(Float_grid_traits()));
  assert(out_segs.size() == 11);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();

  // HPSR Traits
  snap_rounding_2(segs, out_segs, params::output_unique_segments(true).geom_traits(Hot_pixel_traits()));
  snap_rounding_2(segs, out_poly, params::geom_traits(Hot_pixel_traits()));
  assert(out_segs.size() == 8);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();
}

#ifndef CGAL_NO_DEPRECATED_CODE
void test_deprecated(){
  using Deprecated_traits = CGAL::Snap_rounding_traits_2<Kernel>;

  std::vector< Segment_2 > segs;
  std::vector< Polyline_2 > out_poly;
  FT e(std::pow(2, -60));
  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
  segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
  segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));

  // Test new API with deprecated traits
  snap_rounding_2(segs, out_poly, params::geom_traits(Deprecated_traits()));
  assert(out_poly.size() == 7);
  out_poly.clear();

  // Test deprecated API with new traits name
  CGAL::snap_rounding_2<Hot_pixel_traits>(segs.begin(), segs.end(), out_poly, 1.);
  assert(out_poly.size() == 7);
  out_poly.clear();

  // Test deprecated API with deprecated traits
  CGAL::snap_rounding_2<Deprecated_traits>(segs.begin(), segs.end(), out_poly, 1.);
  assert(out_poly.size() == 7);
  out_poly.clear();
}
#endif // CGAL_NO_DEPRECATED_CODE

// Hot-pixel snap rounding of a range of polygons
void test_hot_pixel_polygons() {
  using Polygon_2 = CGAL::Polygon_2<Kernel>;
  std::vector<Polygon_2> in(1), out;
  in[0].push_back(Point_2(0.2, 0.2));
  in[0].push_back(Point_2(10.3, 0.1));
  in[0].push_back(Point_2(10.2, 9.8));
  in[0].push_back(Point_2(0.1, 10.4));

  CGAL::hot_pixel_snap_rounding_2(in, out);
  assert(out.size() == 1);
  const Polygon_2& P = out[0];
  assert(P.size() == 4);
  for (const Point_2& p : {Point_2(0, 0), Point_2(10, 0), Point_2(10, 9), Point_2(0, 10)})
    assert(std::find(P.begin(), P.end(), p) != P.end());
}

// A stateful traits class; it counts the calls to `snap_2_object()`.
struct Counting_traits : public Hot_pixel_traits {
  int* m_calls = nullptr;

  explicit Counting_traits(int* calls) : m_calls(calls) {}

  Snap_2 snap_2_object() const {
    if (m_calls != nullptr) ++*m_calls;
    return Hot_pixel_traits::snap_2_object();
  }
};

// The traits object passed as a named parameter must be the one in use.
void test_traits_instance() {
  std::vector<Segment_2> segs;
  segs.emplace_back(Point_2(0.2, 0.3), Point_2(10.7, 5.1));
  segs.emplace_back(Point_2(0.4, 6.2), Point_2(9.6, 0.8));

  int calls = 0;
  std::vector<Polyline_2> out_poly;
  CGAL::hot_pixel_snap_rounding_2(segs, out_poly, params::geom_traits(Counting_traits(&calls)));
  assert(out_poly.size() == 2);
  assert(calls > 0);

  using Polygon_2 = CGAL::Polygon_2<Kernel>;
  std::vector<Polygon_2> in(1), out;
  for (const Point_2& p : {Point_2(0.2, 0.2), Point_2(10.3, 0.1), Point_2(5.2, 9.8)}) in[0].push_back(p);
  calls = 0;
  CGAL::hot_pixel_snap_rounding_2(in, out, params::geom_traits(Counting_traits(&calls)));
  assert(out.size() == 1);
  assert(calls > 0);
}

int main(int /*argc*/,char** /*argv*/)
{
  test_api();
  test_hot_pixel_polygons();
  test_traits_instance();
#ifndef CGAL_NO_DEPRECATED_CODE
  test_deprecated();
#endif // CGAL_NO_DEPRECATED_CODE

  return(0);
}

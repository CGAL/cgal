#define CGAL_NO_DEPRECATED_WARNING

#include <CGAL/Hot_pixel_snap_rounding_traits_2.h>
#include <CGAL/Double_grid_snap_rounding_traits_2.h>
#include <CGAL/Float_grid_snap_rounding_traits_2.h>
#include <CGAL/Integer_grid_snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <cstdlib>

#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/Snap_rounding_traits_2.h>
#endif // CGAL_NO_DEPRECATED_CODE

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Segment_2                                       Segment_2;
typedef Kernel::Point_2                                         Point_2;
typedef std::vector< Point_2 >                                  Polyline_2;
typedef Kernel::Vector_2                                        Vector_2;
typedef Kernel::FT                                              FT;

void test_api(){
  std::vector< Segment_2 > segs;
  std::vector< Segment_2 > out_segs;
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
  snap_rounding_2(segs, std::back_inserter(out_segs));
  snap_rounding_2(segs, std::back_inserter(out_poly));
  assert(out_segs.size() == 11);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();

  // VSSR Traits
  snap_rounding_2(segs, std::back_inserter(out_segs), CGAL::parameters::geom_traits(CGAL::Float_grid_snap_rounding_traits_2<Kernel>()));
  snap_rounding_2(segs, std::back_inserter(out_poly), CGAL::parameters::geom_traits(CGAL::Float_grid_snap_rounding_traits_2<Kernel>()));
  assert(out_segs.size() == 11);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();

  // HPSR Traits
  snap_rounding_2(segs, std::back_inserter(out_segs), CGAL::parameters::geom_traits(CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>()));
  snap_rounding_2(segs, std::back_inserter(out_poly), CGAL::parameters::geom_traits(CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>()));
  assert(out_segs.size() == 8);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();
}

#ifndef CGAL_NO_DEPRECATED_CODE
void test_deprecated(){
  using Deprecated_traits = CGAL::Snap_rounding_traits_2<Kernel>;
  using New_traits = CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>;

  std::vector< Segment_2 > segs;
  std::vector< Segment_2 > out_segs;
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
  snap_rounding_2(segs, std::back_inserter(out_segs), CGAL::parameters::geom_traits(Deprecated_traits()));
  snap_rounding_2(segs, std::back_inserter(out_poly), CGAL::parameters::geom_traits(Deprecated_traits()));
  assert(out_segs.size() == 8);
  assert(out_poly.size() == 7);
  out_segs.clear();
  out_poly.clear();

  // Test deprecated API with new traits name
  CGAL::snap_rounding_2<New_traits> (segs.begin(), segs.end(), out_poly, 1.);
  assert(out_poly.size() == 7);
  out_poly.clear();

  // Test deprecated API with deprecated traits
  CGAL::snap_rounding_2<Deprecated_traits> (segs.begin(), segs.end(), out_poly, 1.);
  assert(out_poly.size() == 7);
  out_poly.clear();
}
#endif // CGAL_NO_DEPRECATED_CODE

int main(int /*argc*/,char** /*argv*/)
{
  test_api();
#ifndef CGAL_NO_DEPRECATED_CODE
  test_deprecated();
#endif // CGAL_NO_DEPRECATED_CODE

  return(0);
}

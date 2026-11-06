// Test the compile-time detection of geometry traits on curved surfaces used by
// the drawing of arrangements. The drawing code for arrangements on the sphere
// differs from the code for arrangements on the plane. The geometry traits may
// be a (chain of) decorator(s), e.g., Arr_tracing_traits_2 or
// Arr_counting_traits_2, in which case the detection must apply to the last
// (innermost) traits in the chain. Similarly, the geometry traits may be a
// polycurve traits, e.g., Arr_polyline_traits_2 or Arr_polycurve_traits_2,
// whose subcurve traits handle geodesic arcs on the sphere. Decorators and
// polycurve traits can be chained in any order.
// Also, test the construction of the identification curve, which the drawing
// of arrangements on the sphere inserts, through such chains.

#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_counting_traits_2.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_coordinate_converter.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Geodesic_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

// Polycurve traits
using Segment_polyline_traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Segment_polycurve_traits = CGAL::Arr_polycurve_traits_2<Segment_traits>;
using Geodesic_polyline_traits = CGAL::Arr_polyline_traits_2<Geodesic_traits>;
using Geodesic_polycurve_traits = CGAL::Arr_polycurve_traits_2<Geodesic_traits>;
using Nested_geodesic_polycurve_traits =
  CGAL::Arr_polycurve_traits_2<CGAL::Arr_polycurve_traits_2<CGAL::Arr_polycurve_traits_2<Geodesic_traits>>>;

// A traits class derived from the geodesic-arc-on-sphere traits.
class Derived_geodesic_traits : public Geodesic_traits {};

// Geodesic traits with a non-default identification curve
using Rotated_geodesic_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -11, 7>;
using Decorated_rotated_geodesic_traits =
  CGAL::Arr_counting_traits_2<CGAL::Arr_tracing_traits_2<Rotated_geodesic_traits>>;

// The template arguments must be obtained from the last traits in the chain.
static_assert(CGAL::draw_aos::tmpl_args<Decorated_rotated_geodesic_traits>::atan_x == -11);
static_assert(CGAL::draw_aos::tmpl_args<Decorated_rotated_geodesic_traits>::atan_y == 7);
static_assert(std::is_same_v<CGAL::draw_aos::tmpl_args<Decorated_rotated_geodesic_traits>::Kernel, Kernel>);
static_assert(CGAL::draw_aos::tmpl_args<Derived_geodesic_traits>::atan_x == -1);
using Mixed_rotated_geodesic_traits =
  CGAL::Arr_counting_traits_2<CGAL::Arr_polyline_traits_2<CGAL::Arr_tracing_traits_2<
    CGAL::Arr_polycurve_traits_2<Rotated_geodesic_traits>>>>;
static_assert(std::is_same_v<CGAL::draw_aos::innermost_traits_t<Mixed_rotated_geodesic_traits>,
                             Rotated_geodesic_traits>);
static_assert(CGAL::draw_aos::tmpl_args<Mixed_rotated_geodesic_traits>::atan_x == -11);
static_assert(CGAL::draw_aos::tmpl_args<Mixed_rotated_geodesic_traits>::atan_y == 7);

/* Test that the conversion of a point on the sphere to the parameter space
 * with decorated traits is identical to the conversion with the plain traits.
 */
template <typename Decorated, typename Plain>
bool test_coordinate_converter(const std::string& name) {
  Decorated decorated;
  Plain plain;
  CGAL::draw_aos::Arr_coordinate_converter<Decorated> decorated_conv(decorated);
  CGAL::draw_aos::Arr_coordinate_converter<Plain> plain_conv(plain);
  auto p = plain.construct_point_2_object()(0.48, 0.6, 0.64);  // a unit vector
  auto approx_p = plain.approximate_2_object()(p);
  auto decorated_uv = decorated_conv.to_uv(approx_p);
  auto plain_uv = plain_conv.to_uv(approx_p);
  if (decorated_uv == plain_uv) return true;
  std::cerr << "Error: " << name << ": to_uv() = " << decorated_uv << ", expected " << plain_uv << std::endl;
  return false;
}

/* Test that the identification curve constructed through a chain of traits
 * lies on the identification curve of that traits, is vertical, and extends
 * from the south pole (bottom boundary) to the north pole (top boundary).
 */
template <typename Gt>
bool test_identification_curve(const Gt& traits, const std::string& name) {
  auto xcv = CGAL::draw_aos::construct_identification_curve(traits);
  bool on_id = traits.is_on_y_identification_2_object()(xcv);
  bool vertical = traits.is_vertical_2_object()(xcv);
  auto ps_y = traits.parameter_space_in_y_2_object();
  bool ends = (ps_y(xcv, CGAL::ARR_MIN_END) == CGAL::ARR_BOTTOM_BOUNDARY) &&
    (ps_y(xcv, CGAL::ARR_MAX_END) == CGAL::ARR_TOP_BOUNDARY);
  if (on_id && vertical && ends) return true;
  std::cerr << "Error: " << name << ": identification curve: on identification = " << on_id
            << ", vertical = " << vertical << ", from south pole to north pole = " << ends << std::endl;
  return false;
}

template <typename Gt>
bool test(const std::string& name, bool expected) {
  bool agas = CGAL::draw_aos::is_or_derived_from_agas_v<Gt>;
  bool curved = CGAL::draw_aos::is_or_derived_from_curved_surf_traits_v<Gt>;
  if ((agas == expected) && (curved == expected)) return true;
  std::cerr << "Error: " << name << ": is_or_derived_from_agas_v = " << agas
            << ", is_or_derived_from_curved_surf_traits_v = " << curved
            << ", expected " << expected << std::endl;
  return false;
}

int main() {
  bool ok = true;

  // Planar traits
  ok &= test<Segment_traits>("segment traits", false);
  ok &= test<CGAL::Arr_counting_traits_2<Segment_traits>>("counting segment traits", false);
  ok &= test<CGAL::Arr_tracing_traits_2<Segment_traits>>("tracing segment traits", false);
  ok &= test<CGAL::Arr_counting_traits_2<CGAL::Arr_tracing_traits_2<Segment_traits>>>
    ("counting tracing segment traits", false);
  ok &= test<Segment_polyline_traits>("segment polyline traits", false);
  ok &= test<Segment_polycurve_traits>("segment polycurve traits", false);
  ok &= test<CGAL::Arr_counting_traits_2<Segment_polyline_traits>>("counting segment polyline traits", false);

  // Traits on the sphere
  ok &= test<Geodesic_traits>("geodesic traits", true);
  ok &= test<Derived_geodesic_traits>("derived geodesic traits", true);
  ok &= test<CGAL::Arr_counting_traits_2<Geodesic_traits>>("counting geodesic traits", true);
  ok &= test<CGAL::Arr_tracing_traits_2<Geodesic_traits>>("tracing geodesic traits", true);
  ok &= test<CGAL::Arr_counting_traits_2<CGAL::Arr_tracing_traits_2<Geodesic_traits>>>
    ("counting tracing geodesic traits", true);
  ok &= test<CGAL::Arr_tracing_traits_2<CGAL::Arr_counting_traits_2<Derived_geodesic_traits>>>
    ("tracing counting derived geodesic traits", true);
  ok &= test<Geodesic_polyline_traits>("geodesic polyline traits", true);
  ok &= test<Geodesic_polycurve_traits>("geodesic polycurve traits", true);
  ok &= test<CGAL::Arr_counting_traits_2<Geodesic_polyline_traits>>("counting geodesic polyline traits", true);
  ok &= test<CGAL::Arr_polycurve_traits_2<CGAL::Arr_counting_traits_2<Geodesic_traits>>>
    ("polycurve of counting geodesic traits", true);
  ok &= test<Nested_geodesic_polycurve_traits>("nested geodesic polycurve traits", true);
  ok &= test<Mixed_rotated_geodesic_traits>("counting polyline tracing polycurve rotated geodesic traits", true);

  ok &= test_coordinate_converter<CGAL::Arr_counting_traits_2<Geodesic_traits>, Geodesic_traits>
    ("counting geodesic traits");
  ok &= test_coordinate_converter<Derived_geodesic_traits, Geodesic_traits>("derived geodesic traits");
  ok &= test_coordinate_converter<Decorated_rotated_geodesic_traits, Rotated_geodesic_traits>
    ("counting tracing rotated geodesic traits");
  ok &= test_coordinate_converter<Geodesic_polyline_traits, Geodesic_traits>("geodesic polyline traits");
  ok &= test_coordinate_converter<Nested_geodesic_polycurve_traits, Geodesic_traits>
    ("nested geodesic polycurve traits");
  ok &= test_coordinate_converter<Mixed_rotated_geodesic_traits, Rotated_geodesic_traits>
    ("counting polyline tracing polycurve rotated geodesic traits");

  ok &= test_identification_curve(Geodesic_traits(), "geodesic traits");
  ok &= test_identification_curve(Rotated_geodesic_traits(), "rotated geodesic traits");
  ok &= test_identification_curve(CGAL::Arr_counting_traits_2<Geodesic_traits>(), "counting geodesic traits");
  ok &= test_identification_curve(Geodesic_polyline_traits(), "geodesic polyline traits");
  ok &= test_identification_curve(Nested_geodesic_polycurve_traits(), "nested geodesic polycurve traits");
  {
    // Construct the chain from a tracing traits object whose traces are disabled.
    using Tracing_traits = CGAL::Arr_tracing_traits_2<CGAL::Arr_polycurve_traits_2<Rotated_geodesic_traits>>;
    Tracing_traits tracing_traits;
    tracing_traits.disable_all_traces();
    Mixed_rotated_geodesic_traits mixed_traits(&tracing_traits);
    ok &= test_identification_curve(mixed_traits, "counting polyline tracing polycurve rotated geodesic traits");
  }

  if (! ok) return EXIT_FAILURE;
  std::cout << "Passed" << std::endl;
  return EXIT_SUCCESS;
}

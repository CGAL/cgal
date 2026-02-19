#include <vector>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/compute_intersection_polylines.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Random.h>

/*! Convert HSV to RGB color space
 * Converts a given set of HSV values `h', `s', `v' into RGB coordinates.
 * The output RGB values are in the range [0, 255], and the input HSV values
 * are in the ranges h = [0, 360], and s, v = [0, 1], respectively.
 *
 * \param hue Hue component range: [0, 360]
 * \param sat Saturation component range: [0, 1]
 * \param value Value component range: [0, 1]
 * \return tuple<red, green, blue>, where each component is in the range [0, 255]
 */
std::tuple<unsigned char, unsigned char, unsigned char> hsv_to_rgb(double h, double s, double v) {
  double c = v * s;
  double x = c * (1 - std::abs(std::fmod(h / 60.0, 2) - 1));
  double m = v - c;
  double rs, gs, bs;

  if (h < 60)       { rs = c; gs = x; bs = 0; }
  else if (h < 120) { rs = x; gs = c; bs = 0; }
  else if (h < 180) { rs = 0; gs = c; bs = x; }
  else if (h < 240) { rs = 0; gs = x; bs = c; }
  else if (h < 300) { rs = x; gs = 0; bs = c; }
  else              { rs = c; gs = 0; bs = x; }

  unsigned char redc = static_cast<unsigned char>((rs + m) * 255);
  unsigned char greenc = static_cast<unsigned char>((gs + m) * 255);
  unsigned char bluec = static_cast<unsigned char>((bs + m) * 255);
  return std::make_tuple(redc, greenc, bluec);
}

/*!
 */
int main(int argc, char* argv[]) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Point = Kernel::Point_2;
  using Segment = Kernel::Segment_2;
  using Base_traits = CGAL::Arr_segment_traits_2<Kernel>;
  using Traits = CGAL::Arr_curve_data_traits_2<Base_traits, std::size_t>;
  using X_monotone_curve_2 = Traits::X_monotone_curve_2;
  using Polyline = std::vector<std::size_t>;

  std::vector<X_monotone_curve_2> segments = {
    X_monotone_curve_2(Segment(Point(1, 5), Point(8, 5)), 0),
    X_monotone_curve_2(Segment(Point(1, 1), Point(8, 8)), 1),
    X_monotone_curve_2(Segment(Point(3, 1), Point(3, 8)), 2),
    X_monotone_curve_2(Segment(Point(8, 5), Point(8, 8)), 3)
  };

  std::vector<Point> points;
  std::vector<Polyline> polylines;
  Traits traits;
  CGAL::Surface_sweep_2::compute_intersection_polylines(segments.begin(), segments.end(), points, polylines, traits);
  std::size_t i = 0;
  for (auto& polyline : polylines) {
    std::cout << "Polyline " << i++ << ": ";
    for (auto id : polyline) std::cout << id << " ";
    std::cout << std::endl;
  }

  using Arrangement = CGAL::Arrangement_2<Traits>;
  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
  using Edge_const_handle = Arrangement::Halfedge_const_handle;
  using Face_const_handle = Arrangement::Face_const_handle;

  Arrangement arr(&traits);
  CGAL::insert(arr, segments.begin(), segments.end());
  CGAL::Graphics_scene scene;
  CGAL::Graphics_scene_options<Arrangement, Vertex_const_handle, Halfedge_const_handle, Face_const_handle> gso;

  gso.draw_face = [](const Arrangement&, Face_const_handle) { return false; };

  gso.colored_edge = [](const Arrangement&, Edge_const_handle) -> bool { return true; };
  gso.edge_color = [&](const Arrangement& arr, Edge_const_handle eh) -> CGAL::IO::Color {
    auto id = eh->curve().data();
    double h = (static_cast<double>(id) * 360.0) / segments.size();
    double s = 0.8;
    double v = 0.95;
    auto [r, g, b] = hsv_to_rgb(h, s, v);
    return CGAL::IO::Color(r, g, b);
  };

  gso.draw_vertex = [](const Arrangement&, Vertex_const_handle) { return true; }; //! \todo ignored?
  gso.colored_vertex = [](const Arrangement&, Vertex_const_handle) -> bool { return true; };
  gso.vertex_color = [](const Arrangement&, Vertex_const_handle vh) -> CGAL::IO::Color { return CGAL::IO::black(); };

  CGAL::draw(arr, gso, "Arrangement");
  return 0;
}

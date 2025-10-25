#include <string>

#include <CGAL/Random.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/draw_arrangement_2.h>

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
std::tuple<unsigned char, unsigned char, unsigned char> hsv_to_rgb(float hue, float sat, float value) {
  float red, green, blue;
  float fc = value * sat; // Chroma
  float hue_prime = fmod(hue / 60.0f, 6.f);
  float fx = fc * (1.0f - fabs(fmod(hue_prime, 2.f) - 1.f));
  float fm = value - fc;

  if(0 <= hue_prime && hue_prime < 1) {
    red = fc;
    green = fx;
    blue = 0;
  }
  else if(1 <= hue_prime && hue_prime < 2) {
    red = fx;
    green = fc;
    blue = 0;
  }
  else if(2 <= hue_prime && hue_prime < 3) {
    red = 0;
    green = fc;
    blue = fx;
  }
  else if(3 <= hue_prime && hue_prime < 4) {
    red = 0;
    green = fx;
    blue = fc;
  }
  else if(4 <= hue_prime && hue_prime < 5) {
    red = fx;
    green = 0;
    blue = fc;
  }
  else if(5 <= hue_prime && hue_prime < 6) {
    red = fc;
    green = 0;
    blue = fx;
  }
  else {
    red = 0;
    green = 0;
    blue = 0;
  }

  red += fm;
  green += fm;
  blue += fm;

  red *= 255;
  green *= 255;
  blue *= 255;
  unsigned char redc = (unsigned char)red;
  unsigned char greenc = (unsigned char)green;
  unsigned char bluec = (unsigned char)blue;
  return std::make_tuple(redc, greenc, bluec);
}

void draw_rect() {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Arr_non_caching_segment_traits_2<Kernel>;
  using Point = Traits::Point_2;
  using Arrangement = CGAL::Arrangement_2<Traits>;

  Traits traits;
  Arrangement arr(&traits);
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  CGAL::insert(arr, ctr_xcv(Point(-2, -2), Point(2, -2)));
  CGAL::insert(arr, ctr_xcv(Point(2, -2), Point(2, 2)));
  CGAL::insert(arr, ctr_xcv(Point(2, 2), Point(-2, 2)));
  CGAL::insert(arr, ctr_xcv(Point(-2, 2), Point(-2, -2)));

  CGAL::insert(arr, ctr_xcv(Point(-1, -1), Point(1, -1)));
  CGAL::insert(arr, ctr_xcv(Point(1, -1), Point(1, 1)));
  CGAL::insert(arr, ctr_xcv(Point(1, 1), Point(-1, 1)));
  CGAL::insert(arr, ctr_xcv(Point(-1, 1), Point(-1, -1)));

  CGAL::insert(arr, ctr_xcv(Point(-2, -2), Point(-2, -4)));
  CGAL::insert(arr, ctr_xcv(Point(2, -2), Point(4, -2)));

  CGAL::insert(arr, ctr_xcv(Point(0, 0), Point(0, -3)));

  std::cout << arr.number_of_vertices() << ", " << arr.number_of_edges() << ", " << arr.number_of_faces() << std::endl;

  CGAL::Graphics_scene_options<Arrangement, typename Arrangement::Vertex_const_handle,
                               typename Arrangement::Halfedge_const_handle, typename Arrangement::Face_const_handle>
      gso;
  gso.colored_face = [](const Arrangement&, Arrangement::Face_const_handle) -> bool { return true; };

  gso.face_color = [](const Arrangement& arr, Arrangement::Face_const_handle fh) -> CGAL::IO::Color {
    CGAL::Random random((size_t(fh.ptr())));
    float h = 360.0f * random.get_double(0, 1);
    float s = 0.5;
    float v = 0.5;
    auto [r, g, b] = hsv_to_rgb(h, s, v);
    return CGAL::IO::Color(r, g, b);
  };

  CGAL::draw(arr, gso, "rect with hsv colors");
}

void draw_nested() {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Arr_segment_traits_2<Kernel>;
  using Point = Traits::Point_2;
  using Arrangement = CGAL::Arrangement_2<Traits>;
  using X_monotone_curve = Traits::X_monotone_curve_2;

  Arrangement arr;
  auto traits = arr.traits();
  auto ctr_xcv = traits->construct_x_monotone_curve_2_object();

  std::vector<X_monotone_curve> curves;
  {
    // a hexagon centered at the origin
    const double r = 10.0;
    for (int i = 0; i < 6; ++i) {
      int next = (i + 1) % 6;
      Point source(r * cos(i * CGAL_PI / 3), r * sin(i * CGAL_PI / 3));
      Point target(r * cos(next * CGAL_PI / 3), r * sin(next * CGAL_PI / 3));
      curves.push_back(ctr_xcv(source, target));
    }
  }
  {
    // a square inside the hexagon
    const double size = 4.0;
    Point p1(-size, size), p2(size, size), p3(size, -size), p4(-size, -size);
    auto rect = {ctr_xcv(p1, p2), ctr_xcv(p2, p3), ctr_xcv(p3, p4), ctr_xcv(p4, p1)};
    curves.insert(curves.end(), rect.begin(), rect.end());
  }
  {
    // two adjacent triangle inside the square
    Point p1(-1, 0), p2(1, 0), p3(0, sqrt(3)), p4(0, -sqrt(3));
    auto tri1 = {ctr_xcv(p1, p2), ctr_xcv(p2, p3), ctr_xcv(p3, p1)};
    auto tri2 = {ctr_xcv(p1, p2), ctr_xcv(p2, p4), ctr_xcv(p4, p1)};
    curves.insert(curves.end(), tri1.begin(), tri1.end());
    curves.insert(curves.end(), tri2.begin(), tri2.end());
  }
  // a degenerate hole inside the square
  auto degen_seg = ctr_xcv({-1, -3}, {1, -3});
  curves.push_back(degen_seg);
  // an isolated vertex inside the square
  auto iso_point = Point{1, -1};

  CGAL::insert(arr, curves.begin(), curves.end());
  CGAL::insert_point(arr, iso_point);
  CGAL::draw(arr, "nested polygons");
}

void draw_unbounded_linear_grid() {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Arr_linear_traits_2<Kernel>;
  using Point = Traits::Point_2;
  using Segment = Traits::Segment_2;
  using Line = Traits::Line_2;
  using X_monotone_curve = Traits::X_monotone_curve_2;
  using Arrangement = CGAL::Arrangement_2<Traits>;

  Arrangement arr;

  // Insert a n*n grid
  int n = 5;
  for (int i = 0; i < n; ++i) {
    Point p1(i * 5, 0);
    Point p2(i * 5, 1);
    CGAL::insert(arr, X_monotone_curve(Line(p1, p2)));
  }
  for (int i = 0; i < n; ++i) {
    Point p1(0, i * 5);
    Point p2(1, i * 5);
    CGAL::insert(arr, X_monotone_curve(Line(p1, p2)));
  }
  // Generate a inner square(2*2) for all cells
  // And an inner triangle for each square
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      Point p1(i * 5 + 1, j * 5 + 1);
      Point p2(i * 5 + 4, j * 5 + 4);
      CGAL::insert(arr, X_monotone_curve(Segment(p1, Point(p2.x(), p1.y()))));
      CGAL::insert(arr, X_monotone_curve(Segment(Point(p1.x(), p2.y()), p2)));
      CGAL::insert(arr, X_monotone_curve(Segment(p1, Point(p1.x(), p2.y()))));
      CGAL::insert(arr, X_monotone_curve(Segment(Point(p2.x(), p1.y()), p2)));
      // Insert a triangle inside the square
      Point tri_p1(i * 5 + 2, j * 5 + 2);
      Point tri_p2(i * 5 + 3, j * 5 + 2);
      Point tri_p3(i * 5 + 2.5, j * 5 + 3);
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p1, tri_p2)));
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p2, tri_p3)));
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p3, tri_p1)));
      // Connect the triangle to the square
      Point top(i * 5 + 2.5, j * 5 + 4);
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p1, top)));
    }
  }

  CGAL::draw(arr, "unbounded linear grid");
}

void draw_random_segments(int n) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Arr_segment_traits_2<Kernel>;
  using Point = Traits::Point_2;
  using Arrangement = CGAL::Arrangement_2<Traits>;
  using X_monotone_curve = Traits::X_monotone_curve_2;

  Arrangement arr;
  auto traits = arr.traits();
  auto ctr_xcv = traits->construct_x_monotone_curve_2_object();
  CGAL::Random random;

  std::vector<X_monotone_curve> curves;
  for (int i = 0; i < n; ++i) {
    double x1 = random.get_double(-100, 100);
    double y1 = random.get_double(-100, 100);
    double x2 = random.get_double(-100, 100);
    double y2 = random.get_double(-100, 100);
    curves.push_back(ctr_xcv(Point(x1, y1), Point(x2, y2)));
  }
  CGAL::insert(arr, curves.begin(), curves.end());
  CGAL::draw(arr, (std::to_string(n) + " random segments").c_str());
}

int main() {
  draw_rect();
  draw_nested();
  draw_unbounded_linear_grid();
  draw_random_segments(100);
  return EXIT_SUCCESS;
}

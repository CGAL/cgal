#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/draw_arrangement_2.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<Kernel>;
using Point = Traits::Point_2;
using Arrangement_2 = CGAL::Arrangement_2<Traits>;

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
std::tuple<float, float, float>
hsv_to_rgb(float hue, float sat, float value) {
  float red, green, blue;
  float fc = value * sat; // Chroma
  float hue_prime = fmod(hue / 60.0, 6);
  float fx = fc * (1.0 - fabs(fmod(hue_prime, 2) - 1.0));
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
  return std::make_tuple(red, green, blue);
}

int main() {
  Traits traits;
  Arrangement_2 arr(&traits);
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  CGAL::insert(arr, ctr_xcv(Point(-2,-2), Point(2,-2)));
  CGAL::insert(arr, ctr_xcv(Point(2,-2), Point(2,2)));
  CGAL::insert(arr, ctr_xcv(Point(2,2), Point(-2,2)));
  CGAL::insert(arr, ctr_xcv(Point(-2,2), Point(-2,-2)));

  CGAL::insert(arr, ctr_xcv(Point(-1,-1), Point(1,-1)));
  CGAL::insert(arr, ctr_xcv(Point(1,-1), Point(1,1)));
  CGAL::insert(arr, ctr_xcv(Point(1,1), Point(-1,1)));
  CGAL::insert(arr, ctr_xcv(Point(-1,1), Point(-1,-1)));

  CGAL::insert(arr, ctr_xcv(Point(-2,-2), Point(-2,-4)));
  CGAL::insert(arr, ctr_xcv(Point(2,-2), Point(4,-2)));

  CGAL::insert(arr, ctr_xcv(Point(0,0), Point(0,-3)));

  std::cout << arr.number_of_vertices() << ", "
            << arr.number_of_edges() << ", "
            << arr.number_of_faces() << std::endl;

  std::size_t id(0);
  CGAL::draw(arr, [&] (Arrangement_2::Face_const_handle) -> CGAL::IO::Color {
                    float h = 360.0 * id++ / arr.number_of_faces();
                    float s = 0.5;
                    float v = 0.5;
                    float r, g, b;
                    std::tie(r, g, b) = hsv_to_rgb(h, s, v);
                    return CGAL::IO::Color(r, g, b);
                  }, "hsv colors", true);

  return EXIT_SUCCESS;
}

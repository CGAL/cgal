#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/draw_arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
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
std::tuple<unsigned char, unsigned char, unsigned char>
hsv_to_rgb(float hue, float sat, float value) {
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

 CGAL::Graphics_scene_options<Arrangement_2,
                              typename Arrangement_2::Vertex_const_handle,
                              typename Arrangement_2::Halfedge_const_handle,
                              typename Arrangement_2::Face_const_handle> gso;
 gso.colored_face=[](const Arrangement_2&, Arrangement_2::Face_const_handle) -> bool
 { return true; };

 gso.face_color=[&id](const Arrangement_2& arr, Arrangement_2::Face_const_handle) -> CGAL::IO::Color
 {
   float h = 360.0f * id++ / arr.number_of_faces();
   float s = 0.5;
   float v = 0.5;
   auto [r, g, b] = hsv_to_rgb(h, s, v);
   return CGAL::IO::Color(r,g,b);
 };

 CGAL::draw(arr, gso, "hsv colors");

 return EXIT_SUCCESS;
}

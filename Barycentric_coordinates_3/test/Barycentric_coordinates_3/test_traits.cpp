#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>
#include "include/utils.h"

struct Custom_barycentric_traits
{
  // New requirements {
  struct FT
  {
    FT() {}
    FT(double) {}
    FT operator/(FT)const { return FT(); }
    FT operator/=(FT) { return FT(); }
    FT operator*=(FT) { return FT(); }
    FT operator-(const FT)const { return FT(); }
    FT operator+(FT)const { return FT(); }
    FT operator*(FT)const { return FT(); }
    bool operator<=(FT)const { return false; }
    bool operator==(FT)const { return false; }
    bool operator>=(FT)const { return false; }
    bool operator!=(FT)const { return false; }
    bool operator<(FT)const { return false; }
    bool operator>(FT)const { return false; }
    FT& operator+=(FT) { return *this; }
  };

  struct Point_3
  {
    Point_3() {}
    Point_3(FT, FT, FT) {}
  };

  struct Point_2
  {
    Point_2() {}
    Point_2(FT, FT) {}
  };

  struct Vector_3 {
  };

  friend Vector_3 operator*(const int, const Vector_3) {return Vector_3();}

  struct Vector_2 {
  };

  struct Triangle_2 {
  };

  struct Iso_rectangle_2 {
  };

  struct Plane_3 {
    Plane_3() {};
  };

  struct Compare_x_2 {
    CGAL::Comparison_result operator()(const Point_2& p, const Point_2& q) { return CGAL::EQUAL; }
  };

  struct Compare_y_2 {
    CGAL::Comparison_result operator()(const Point_2& p, const Point_2& q) { return CGAL::EQUAL; }
  };

  struct Compute_approximate_angle_3 {
    FT operator()(const Vector_3&, const Vector_3&) { return FT(); }
  };

  struct Compute_area_2 {
    FT operator()(const Point_2& p, const Point_2& q, const Point_2& r) { return FT(); }
    FT operator()(const Iso_rectangle_2& r) { return FT(); }
    FT operator()(const Triangle_2& t) { return FT(); }
  };

  struct Compute_scalar_product_3 {
    FT operator()(const Vector_3&, const Vector_3&) { return FT(); }
  };

  struct Construct_cross_product_vector_3 {
    Vector_3 operator()(const Vector_3&, const Vector_3&) { return Vector_3(); }
  };

  struct Construct_divided_vector_3 {
    Vector_3 operator()(const Vector_3&, const FT&) { return Vector_3(); }
  };

  struct Construct_plane_3 {
    Plane_3 operator()(const Point_3&, const Point_3&, const Point_3&) { return Plane_3(); }
  };

  struct Construct_projected_xy_point_2 {
    Point_2 operator()(const Plane_3&, const Point_3&) { return Point_2(); }
  };

  struct Compute_determinant_3 {
    FT operator()(const Vector_3&, const Vector_3&, const Vector_3&) { return FT(); }
  };

  struct Compute_scalar_product_2 {
    FT operator()(const Vector_2&, const Vector_2&) { return FT(); }
  };

  struct Compute_squared_distance_2 {
    FT operator()(const Point_2&, const Point_2&) { return FT(); }
  };

  struct Compute_squared_length_3 {
    FT operator()(const Vector_3&) { return FT(); }
  };

  struct Compute_volume_3 {
    FT operator()(const Point_3&, const Point_3&, const Point_3&, const Point_3&) const { return FT(); }
  };

  struct Construct_vector_2 {
    Vector_2         operator()(const Point_2&, const Point_2&) const { return Vector_2(); }
  };

  struct Construct_vector_3 {
    Vector_3         operator()(const Point_3&, const Point_3&) const { return Vector_3(); }
  };

  struct Collinear_2 {
    bool operator()(const Point_2&, const Point_2&, const Point_2&) const { return true; }
  };

  struct Collinear_are_ordered_along_line_2 {
    bool operator()(const Point_2&, const Point_2&, const Point_2&) const { return true; }
  };

  struct Equal_2 {
    bool operator()(const Point_2& p, const Point_2& q) const { return true; }
  };

  struct Orientation_2 {
    CGAL::Orientation operator()(const Point_2&, const Point_2&, const Point_2&) const { return CGAL::COLLINEAR; }
  };

  struct Less_xy_2 {
    bool operator()(const Point_2& p, const Point_2& q) const { return true; }
  };

  Collinear_2     collinear_2_object() const { return Collinear_2(); }
  Collinear_are_ordered_along_line_2 collinear_are_ordered_along_line_2_object() const { return Collinear_are_ordered_along_line_2(); }
  Compare_x_2       compare_x_2_object() const { return Compare_x_2(); }
  Compare_y_2       compare_y_2_object() const { return Compare_y_2(); }
  Compute_approximate_angle_3 compute_approximate_angle_3_object() const { return Compute_approximate_angle_3(); }
  Compute_area_2      compute_area_2_object() const { return Compute_area_2(); }
  Compute_scalar_product_3    compute_scalar_product_3_object() const { return Compute_scalar_product_3(); }
  Compute_scalar_product_2    compute_scalar_product_2_object() const { return Compute_scalar_product_2(); }
  Compute_squared_length_3    compute_squared_length_3_object() const { return Compute_squared_length_3(); }
  Compute_squared_distance_2  compute_squared_distance_2_object() const { return Compute_squared_distance_2(); }
  Construct_cross_product_vector_3 construct_cross_product_vector_3_object() const { return Construct_cross_product_vector_3(); }
  Compute_determinant_3 compute_determinant_3_object() const { return Compute_determinant_3(); }
  Compute_volume_3    compute_volume_3_object() const { return Compute_volume_3(); }
  Construct_divided_vector_3 construct_divided_vector_3_object() const { return Construct_divided_vector_3(); }
  Construct_plane_3   construct_plane_3_object() const { return Construct_plane_3(); }
  Construct_projected_xy_point_2 construct_projected_xy_point_2_object() const { return Construct_projected_xy_point_2(); }
  Construct_vector_2  construct_vector_2_object() const { return Construct_vector_2(); }
  Construct_vector_3  construct_vector_3_object() const { return Construct_vector_3(); }
  Equal_2            equal_2_object() const { return Equal_2(); }
  Orientation_2      orientation_2_object() const { return Orientation_2(); }
  Less_xy_2          less_xy_2_object() const { return Less_xy_2(); }

  static Custom_barycentric_traits instance() { return Custom_barycentric_traits(); }
private:
  Custom_barycentric_traits() {};
};

namespace CGAL {
Custom_barycentric_traits::FT abs(const Custom_barycentric_traits::FT&) {
  return Custom_barycentric_traits::FT();
}

Custom_barycentric_traits::FT approximate_sqrt(const Custom_barycentric_traits::FT&) {
  return Custom_barycentric_traits::FT();
}

bool is_zero(const Custom_barycentric_traits::FT&) {
  return true;
}
Custom_barycentric_traits::FT sqrt(const Custom_barycentric_traits::FT&) {
  return Custom_barycentric_traits::FT();
}
Custom_barycentric_traits::FT to_double(const Custom_barycentric_traits::FT&) {
  return Custom_barycentric_traits::FT();
}

template<>
class Kernel_traits<Custom_barycentric_traits::Point_3> {
public:
  using Kernel = Custom_barycentric_traits;
};

template<>
class Kernel_traits<Custom_barycentric_traits::Point_2> {
public:
  using Kernel = Custom_barycentric_traits;
};

} // namespace CGAL

// Typedefs.
using Mesh = CGAL::Surface_mesh<Custom_barycentric_traits::Point_3>;
using MV = CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Mesh, Custom_barycentric_traits>;
using WP = CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Custom_barycentric_traits>;
using DH = CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Custom_barycentric_traits>;

int main() {
  Mesh tetrahedron;
  std::vector<Custom_barycentric_traits::Point_3> vertices;
  std::tie(tetrahedron, vertices) = tests::get_irregular_tetrahedron<Custom_barycentric_traits, Mesh>();

  std::vector<Custom_barycentric_traits::FT> coordinates;

  MV mv(tetrahedron, CGAL::Barycentric_coordinates::Computation_policy_3::FAST_WITH_EDGE_CASES, Custom_barycentric_traits::instance());
  WP wp(tetrahedron, CGAL::Barycentric_coordinates::Computation_policy_3::FAST_WITH_EDGE_CASES, Custom_barycentric_traits::instance());
  DH dh(tetrahedron, CGAL::Barycentric_coordinates::Computation_policy_3::FAST_WITH_EDGE_CASES, Custom_barycentric_traits::instance());

  mv(Custom_barycentric_traits::Point_3(), std::back_inserter(coordinates));
  wp(Custom_barycentric_traits::Point_3(), std::back_inserter(coordinates));
  dh(Custom_barycentric_traits::Point_3(), std::back_inserter(coordinates));

  CGAL::Barycentric_coordinates::wachspress_coordinates_3(tetrahedron, Custom_barycentric_traits::Point_3(), std::back_inserter(coordinates), CGAL::parameters::geom_traits(Custom_barycentric_traits::instance()));
  CGAL::Barycentric_coordinates::discrete_harmonic_coordinates_3(tetrahedron, Custom_barycentric_traits::Point_3(), std::back_inserter(coordinates), CGAL::parameters::geom_traits(Custom_barycentric_traits::instance()));
  CGAL::Barycentric_coordinates::mean_value_coordinates_3(tetrahedron, Custom_barycentric_traits::Point_3(), std::back_inserter(coordinates), CGAL::parameters::geom_traits(Custom_barycentric_traits::instance()));


  return EXIT_SUCCESS;
}

#include <CGAL/Arr_counting_traits_2.h>

#include "decorator_test.h"

using namespace CGAL;

template <typename Ct>
class Counting_traits_test {
public:
  using Operation_id = typename Ct::Operation_id;

  Counting_traits_test(const Ct& ct) : m_ct(ct) {}

  template <typename Lambda>
  void test_counting(Operation_id op_id, Lambda fn) {
    std::size_t before = m_ct.count(op_id);
    fn();
    std::size_t after = m_ct.count(op_id);
    if (after == before + 1) return;
    m_ct.print(std::cerr, op_id);
    std::cerr << "Expected the above counter to increase by 1, got: before = " << before << ", after = " << after
              << std::endl;
    exit(1);
  }

  const Ct& m_ct;
};

template <typename Ct>
void test_counting_traits(const Decorating_traits_test_objects<Ct>& objs) {
  using X_monotone_curve_2 = typename Ct::X_monotone_curve_2;
  using Point_2 = typename Ct::Point_2;

  const auto& ct = objs.traits;
  Counting_traits_test<Ct> t(ct);
  const auto& pt = objs.xcv_seg_pt;
  const auto& pt_lr = objs.pt_lr;
  const auto& cv = objs.cv;
  const auto& xcv_seg = objs.xcv_seg;
  auto seg_min = ct.construct_min_vertex_2_object()(xcv_seg);
  auto seg_max = ct.construct_max_vertex_2_object()(xcv_seg);
  const auto& xcv_lr = objs.xcv_lr;
  auto xcv_lr_end = objs.xcv_lr_end;
  const auto& xcv_bt = objs.xcv_bt;
  auto xcv_bt_end = objs.xcv_bt_end;

  if constexpr (has_compare_x_2<Ct>::value)
    t.test_counting(Ct::COMPARE_X_2_OP, [&]() { ct.compare_x_2_object()(pt, pt); });
  if constexpr (has_compare_xy_2<Ct>::value)
    t.test_counting(Ct::COMPARE_XY_2_OP, [&]() { ct.compare_xy_2_object()(pt, pt); });
  if constexpr (has_construct_min_vertex_2<Ct>::value)
    t.test_counting(Ct::CONSTRUCT_MIN_VERTEX_2_OP, [&]() { ct.construct_min_vertex_2_object()(xcv_seg); });
  if constexpr (has_construct_max_vertex_2<Ct>::value)
    t.test_counting(Ct::CONSTRUCT_MAX_VERTEX_2_OP, [&]() { ct.construct_max_vertex_2_object()(xcv_seg); });
  if constexpr (has_is_vertical_2<Ct>::value)
    t.test_counting(Ct::IS_VERTICAL_2_OP, [&]() { ct.is_vertical_2_object()(xcv_seg); });
  if constexpr (has_compare_y_at_x_2<Ct>::value)
    t.test_counting(Ct::COMPARE_Y_AT_X_2_OP, [&]() { ct.compare_y_at_x_2_object()(pt, xcv_seg); });
  if constexpr (has_equal_2<Ct>::value) {
    t.test_counting(Ct::EQUAL_2_POINTS_OP, [&]() { ct.equal_2_object()(pt, pt); });
    t.test_counting(Ct::EQUAL_2_CURVES_OP, [&]() { ct.equal_2_object()(xcv_seg, xcv_seg); });
  }
  if constexpr (has_compare_y_at_x_left_2<Ct>::value){
    t.test_counting(Ct::COMPARE_Y_AT_X_LEFT_2_OP,
                    [&]() { ct.compare_y_at_x_left_2_object()(xcv_seg, xcv_seg, seg_max); });
  }
  if constexpr (has_compare_y_at_x_right_2<Ct>::value) {
    t.test_counting(Ct::COMPARE_Y_AT_X_RIGHT_2_OP,
                    [&]() { ct.compare_y_at_x_right_2_object()(xcv_seg, xcv_seg, seg_min); });
  }
  if constexpr (has_parameter_space_in_x_2<Ct>::value && ! is_or_derived_from_conic_traits_v<unwrap_t<Ct>>) {
    t.test_counting(Ct::PARAMETER_SPACE_IN_X_2_POINT_OP, [&]() { ct.parameter_space_in_x_2_object()(pt); });
    t.test_counting(Ct::PARAMETER_SPACE_IN_X_2_CURVE_END_OP,
                    [&]() { ct.parameter_space_in_x_2_object()(xcv_seg, ARR_MIN_END); });
  }
  if constexpr (has_parameter_space_in_y_2<Ct>::value && ! is_or_derived_from_conic_traits_v<unwrap_t<Ct>>) {
    t.test_counting(Ct::PARAMETER_SPACE_IN_Y_2_POINT_OP, [&]() { ct.parameter_space_in_y_2_object()(pt); });
    t.test_counting(Ct::PARAMETER_SPACE_IN_Y_2_CURVE_END_OP,
                    [&]() { ct.parameter_space_in_y_2_object()(xcv_seg, ARR_MIN_END); });
  }
  if constexpr (has_make_x_monotone_2<Ct>::value) {
    t.test_counting(Ct::MAKE_X_MONOTONE_2_OP, [&]() {
      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;
      std::vector<Make_x_monotone_result> output;
      ct.make_x_monotone_2_object()(cv, std::back_inserter(output));
    });
  }
  X_monotone_curve_2 xcv1, xcv2, xcv_merged; // for split/merge tests
  if constexpr (has_split_2<Ct>::value)
    t.test_counting(Ct::SPLIT_2_OP, [&]() { ct.split_2_object()(xcv_seg, pt, xcv1, xcv2); });
  if constexpr (has_are_mergeable_2<Ct>::value)
    t.test_counting(Ct::ARE_MERGEABLE_2_OP, [&]() { ct.are_mergeable_2_object()(xcv1, xcv2); });
  if constexpr (has_merge_2<Ct>::value)
    t.test_counting(Ct::MERGE_2_OP, [&]() { ct.merge_2_object()(xcv1, xcv2, xcv_merged); });
  if constexpr (has_intersect_2<Ct>::value) {
    t.test_counting(Ct::INTERSECT_2_OP, [&]() {
      using Multiplicity = typename Ct::Multiplicity;
      using Intersection_result = std::variant<std::pair<Point_2, Multiplicity>, X_monotone_curve_2>;
      std::vector<Intersection_result> inters;
      ct.intersect_2_object()(xcv_seg, xcv_seg, std::back_inserter(inters));
    });
  }
  if constexpr (has_construct_opposite_2<Ct>::value)
    t.test_counting(Ct::CONSTRUCT_2_OPPOSITE_2_OP, [&]() { ct.construct_opposite_2_object()(xcv_seg); });
  if constexpr (has_compare_endpoints_xy_2<Ct>::value)
    t.test_counting(Ct::COMPARE_ENDPOINTS_XY_2_OP, [&]() { ct.compare_endpoints_xy_2_object()(xcv_seg); });
  if constexpr (has_approximate_2<Ct>::value) {
    t.test_counting(Ct::APPROXIMATE_2_COORD_OP, [&]() { ct.approximate_2_object()(pt, 0); });
    if constexpr (has_approximate_2_point<Ct>::value)
      t.test_counting(Ct::APPROXIMATE_2_POINT_OP, [&]() { ct.approximate_2_object()(pt); });
    if constexpr (has_approximate_2_xcv<Ct>::value) {
      t.test_counting(Ct::APPROXIMATE_2_CURVE_OP, [&]() {
        using Ct_approximate_2 = typename Ct::Approximate_2;
        using Approximate_point_2 = typename Ct_approximate_2::Approximate_point_2;
        std::vector<Approximate_point_2> approx_points;
        ct.approximate_2_object()(xcv_seg, 0.5, std::back_inserter(approx_points));
      });
    }
    if constexpr (has_approximate_2_xcv_bounds<Ct>::value) {
      t.test_counting(Ct::APPROXIMATE_2_BOUNDED_CURVE_OP, [&]() {
        using Ct_approximate_2 = typename Ct::Approximate_2;
        using Approximate_point_2 = typename Ct_approximate_2::Approximate_point_2;
        std::vector<Approximate_point_2> approx_points;
        ct.approximate_2_object()(xcv_seg, 0.5, std::back_inserter(approx_points), Bbox_2(0, 0, 1, 1));
      });
    }
  }
  if constexpr (has_compare_x_on_boundary_2<Ct>::value) {
    t.test_counting(Ct::COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP,
                    [&]() { ct.compare_x_on_boundary_2_object()(pt, xcv_bt, xcv_bt_end); });
    t.test_counting(Ct::COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP,
                    [&]() { ct.compare_x_on_boundary_2_object()(xcv_bt, xcv_bt_end, xcv_bt, xcv_bt_end); });
  }
  if constexpr (has_compare_x_near_boundary_2<Ct>::value) {
    t.test_counting(Ct::COMPARE_X_NEAR_BOUNDARY_2_OP,
                    [&]() { ct.compare_x_near_boundary_2_object()(xcv_bt, xcv_bt, xcv_bt_end); });
  }
  if constexpr (has_compare_y_on_boundary_2<Ct>::value) {
    t.test_counting(Ct::COMPARE_Y_ON_BOUNDARY_2_OP,
                    [&]() { ct.compare_y_on_boundary_2_object()(pt_lr, pt_lr); });
  }
  if constexpr (has_compare_y_near_boundary_2<Ct>::value) {
    t.test_counting(Ct::COMPARE_Y_NEAR_BOUNDARY_2_OP,
                    [&]() { ct.compare_y_near_boundary_2_object()(xcv_lr, xcv_lr, xcv_lr_end); });
  }
  if constexpr (has_is_on_x_identification_2<Ct>::value) {
    t.test_counting(Ct::IS_ON_X_IDENTIFICATION_2_POINT_OP, [&]() { ct.is_on_x_identification_2_object()(pt); });
    t.test_counting(Ct::IS_ON_X_IDENTIFICATION_2_CURVE_OP, [&]() { ct.is_on_x_identification_2_object()(xcv_seg); });
  }
  if constexpr (has_is_on_y_identification_2<Ct>::value) {
    t.test_counting(Ct::IS_ON_Y_IDENTIFICATION_2_POINT_OP, [&]() { ct.is_on_y_identification_2_object()(pt); });
    t.test_counting(Ct::IS_ON_Y_IDENTIFICATION_2_CURVE_OP, [&]() { ct.is_on_y_identification_2_object()(xcv_seg); });
  }
}

void test_counting_algebraic_traits() {
  using Base_traits = Arr_algebraic_segment_traits_2<CORE::BigInt>;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Counting_traits, Base_traits>();

  Algebraic_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_bezier_curve_traits() {
  using Base_traits = Bezier_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Bezier_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_circle_segment_traits() {
  using Base_traits = Circle_segment_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Circle_segment_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_conic_traits() {
  using Base_traits = Conic_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Conic_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_consolidated_curve_data_traits() {
  using Base_traits = Consolidated_curve_data_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Consolidated_curve_data_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_curve_data_traits() {
  using Base_traits = Curve_data_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Curve_data_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_geodesic_arc_on_sphere_traits() {
  using Base_traits = Arr_geodesic_arc_on_sphere_traits_2<Exact_predicates_exact_constructions_kernel>;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Counting_traits, Base_traits>();

  Geodesic_arc_on_sphere_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_linear_traits() {
  using Base_traits = Linear_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Counting_traits, Base_traits>();

  Linear_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_non_caching_segment_traits() {
  using Base_traits = Non_caching_segment_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Non_caching_segment_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

// void test_counting_polyline_traits() {
//   using Base_traits = Polyline_traits_decorating_test_types::Base_traits;
//   using Counting_traits = Arr_counting_traits_2<Base_traits>;
//   test_type_consistency<Base_traits, Counting_traits>();

//   Polyline_traits_decorating_test_objects<Counting_traits> objs;
//   test_counting_traits(objs);
// }

void test_counting_rational_funcion_traits() {
  using Base_traits = Rational_function_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Rational_function_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

void test_counting_segment_traits() {
  using Base_traits = typename Segment_traits_decorating_test_types::Base_traits;
  using Counting_traits = Arr_counting_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Counting_traits>();

  Segment_traits_decorating_test_objects<Counting_traits> objs;
  test_counting_traits(objs);
}

int main() {
  test_counting_algebraic_traits();
  test_counting_bezier_curve_traits();
  test_counting_circle_segment_traits();
  test_counting_conic_traits();
  test_counting_consolidated_curve_data_traits();
  test_counting_curve_data_traits();
  test_counting_geodesic_arc_on_sphere_traits();
  test_counting_linear_traits();
  test_counting_non_caching_segment_traits();
  // test_counting_polyline_traits();
  test_counting_rational_funcion_traits();
  test_counting_segment_traits();
  return 0;
}

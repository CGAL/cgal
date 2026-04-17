#include <iostream>
#include <sstream>
#include <streambuf>

#include <CGAL/Arr_tracing_traits_2.h>

#include "decorator_test.h"

template <typename Tt>
class Tracing_traits_test {
private:
  static std::streambuf* redirect_cout(std::ostringstream& rd_stream) {
    std::cout.flush();
    std::streambuf* old_buf = std::cout.rdbuf();
    rd_stream.str("");
    std::cout.rdbuf(rd_stream.rdbuf());
    return old_buf;
  }

  static void restore_cout(std::streambuf* old_buf) {
    std::cout.flush();
    std::cout.rdbuf(old_buf);
  }

public:
  using Operation_id = typename Tt::Operation_id;

  Tracing_traits_test(Tt& ct) : m_tt(ct)
  { m_old_buf = redirect_cout(m_rd_stream); }

  ~Tracing_traits_test()
  { restore_cout(m_old_buf); }

  template <typename Lambda>
  void test_tracing(Operation_id op_id, Lambda fn) {
    m_tt.enable_trace(op_id);
    fn();
    if (m_rd_stream.str().empty()) {
      std::cerr << "Expected tracing output for operation id " << op_id << ", got none." << std::endl;
      exit(1);
    }

    m_tt.disable_trace(op_id);
    m_rd_stream.str("");
    fn();
    if (!m_rd_stream.str().empty()) {
      std::cerr << "Did not expect tracing output for operation id " << op_id << ", got some." << std::endl;
      exit(1);
    }
  }

  Tt& m_tt;
  std::ostringstream m_rd_stream;
  std::streambuf* m_old_buf;
};

template <typename Tt>
void test_tracing_traits(Decorating_traits_test_objects<Tt>& objs) {
  using X_monotone_curve_2 = typename Tt::X_monotone_curve_2;
  using Point_2 = typename Tt::Point_2;

  Tt& tt = objs.traits;
  Tracing_traits_test<Tt> t(tt);
  const auto& pt = objs.xcv_seg_pt;
  const auto& pt_lr = objs.pt_lr;
  const auto& pt_bt = objs.pt_bt;
  const auto& cv = objs.cv;
  const auto& xcv_seg = objs.xcv_seg;
  auto seg_min = tt.construct_min_vertex_2_object()(xcv_seg);
  auto seg_max = tt.construct_max_vertex_2_object()(xcv_seg);
  const auto& xcv_lr = objs.xcv_lr;
  auto xcv_lr_end = objs.xcv_lr_end;
  const auto& xcv_bt = objs.xcv_bt;
  auto xcv_bt_end = objs.xcv_bt_end;

  if constexpr (has_compare_x_2<Tt>::value)
    t.test_tracing(Tt::COMPARE_X_2_OP, [&]() { tt.compare_x_2_object()(pt, pt); });
  if constexpr (has_compare_xy_2<Tt>::value)
    t.test_tracing(Tt::COMPARE_XY_2_OP, [&]() { tt.compare_xy_2_object()(pt, pt); });
  if constexpr (has_construct_min_vertex_2<Tt>::value)
    t.test_tracing(Tt::CONSTRUCT_MIN_VERTEX_2_OP, [&]() { tt.construct_min_vertex_2_object()(xcv_seg); });
  if constexpr (has_construct_max_vertex_2<Tt>::value)
    t.test_tracing(Tt::CONSTRUCT_MAX_VERTEX_2_OP, [&]() { tt.construct_max_vertex_2_object()(xcv_seg); });
  if constexpr (has_is_vertical_2<Tt>::value)
    t.test_tracing(Tt::IS_VERTICAL_2_OP, [&]() { tt.is_vertical_2_object()(xcv_seg); });
  if constexpr (has_compare_y_at_x_2<Tt>::value)
    t.test_tracing(Tt::COMPARE_Y_AT_X_2_OP, [&]() { tt.compare_y_at_x_2_object()(pt, xcv_seg); });
  if constexpr (has_equal_2<Tt>::value) {
    t.test_tracing(Tt::EQUAL_2_POINTS_OP, [&]() { tt.equal_2_object()(pt, pt); });
    t.test_tracing(Tt::EQUAL_2_CURVES_OP, [&]() { tt.equal_2_object()(xcv_seg, xcv_seg); });
  }
  if constexpr (has_compare_y_at_x_left_2<Tt>::value){
    t.test_tracing(Tt::COMPARE_Y_AT_X_LEFT_2_OP,
                   [&]() { tt.compare_y_at_x_left_2_object()(xcv_seg, xcv_seg, seg_max); });
  }
  if constexpr (has_compare_y_at_x_right_2<Tt>::value) {
    t.test_tracing(Tt::COMPARE_Y_AT_X_RIGHT_2_OP,
                   [&]() { tt.compare_y_at_x_right_2_object()(xcv_seg, xcv_seg, seg_min); });
  }
  if constexpr (has_parameter_space_in_x_2<Tt>::value && ! is_or_derived_from_conic_traits_v<unwrap_t<Tt>>) {
    t.test_tracing(Tt::PARAMETER_SPACE_IN_X_2_OP, [&]() { tt.parameter_space_in_x_2_object()(pt); });
    t.test_tracing(Tt::PARAMETER_SPACE_IN_X_2_OP,
                   [&]() { tt.parameter_space_in_x_2_object()(xcv_seg, ARR_MIN_END); });
  }
  if constexpr (has_parameter_space_in_y_2<Tt>::value && ! is_or_derived_from_conic_traits_v<unwrap_t<Tt>>) {
    t.test_tracing(Tt::PARAMETER_SPACE_IN_Y_2_OP, [&]() { tt.parameter_space_in_y_2_object()(pt); });
    t.test_tracing(Tt::PARAMETER_SPACE_IN_Y_2_OP,
                   [&]() { tt.parameter_space_in_y_2_object()(xcv_seg, ARR_MIN_END); });
  }
  if constexpr (has_make_x_monotone_2<Tt>::value) {
    t.test_tracing(Tt::MAKE_X_MONOTONE_2_OP, [&]() {
      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;
      std::vector<Make_x_monotone_result> output;
      tt.make_x_monotone_2_object()(cv, std::back_inserter(output));
    });
  }
  X_monotone_curve_2 xcv1, xcv2, xcv_merged; // for split/merge tests
  if constexpr (has_split_2<Tt>::value)
    t.test_tracing(Tt::SPLIT_2_OP, [&]() { tt.split_2_object()(xcv_seg, pt, xcv1, xcv2); });
  if constexpr (has_are_mergeable_2<Tt>::value)
    t.test_tracing(Tt::ARE_MERGEABLE_2_OP, [&]() { tt.are_mergeable_2_object()(xcv1, xcv2); });
  if constexpr (has_merge_2<Tt>::value)
    t.test_tracing(Tt::MERGE_2_OP, [&]() { tt.merge_2_object()(xcv1, xcv2, xcv_merged); });
  if constexpr (has_intersect_2<Tt>::value) {
    t.test_tracing(Tt::INTERSECT_2_OP, [&]() {
      using Multiplicity = typename Tt::Multiplicity;
      using Intersection_result = std::variant<std::pair<Point_2, Multiplicity>, X_monotone_curve_2>;
      std::vector<Intersection_result> inters;
      tt.intersect_2_object()(xcv_seg, xcv_seg, std::back_inserter(inters));
    });
  }
  if constexpr (has_construct_opposite_2<Tt>::value)
    t.test_tracing(Tt::CONSTRUCT_2_OPPOSITE_2_OP, [&]() { tt.construct_opposite_2_object()(xcv_seg); });
  if constexpr (has_compare_endpoints_xy_2<Tt>::value)
    t.test_tracing(Tt::COMPARE_ENDPOINTS_XY_2_OP, [&]() { tt.compare_endpoints_xy_2_object()(xcv_seg); });
  if constexpr (has_approximate_2<Tt>::value) {
    t.test_tracing(Tt::APPROXIMATE_2_OP, [&]() { tt.approximate_2_object()(pt, 0); });
    if constexpr (has_approximate_2_point<Tt>::value)
      t.test_tracing(Tt::APPROXIMATE_2_OP, [&]() { tt.approximate_2_object()(pt); });
    if constexpr (has_approximate_2_xcv<Tt>::value) {
      t.test_tracing(Tt::APPROXIMATE_2_OP, [&]() {
        using Approximate_point_2 = typename Tt::Approximate_point_2;
        std::vector<Approximate_point_2> approx_points;
        tt.approximate_2_object()(xcv_seg, 0.5, std::back_inserter(approx_points));
      });
    }
    if constexpr (has_approximate_2_xcv_bounds<Tt>::value) {
      t.test_tracing(Tt::APPROXIMATE_2_OP, [&]() {
        using Approximate_point_2 = typename Tt::Approximate_point_2;
        std::vector<Approximate_point_2> approx_points;
        tt.approximate_2_object()(xcv_seg, 0.5, std::back_inserter(approx_points), Bbox_2(0, 0, 1, 1));
      });
    }
  }
  if constexpr (has_compare_x_on_boundary_2<Tt>::value) {
    t.test_tracing(Tt::COMPARE_X_ON_BOUNDARY_2_OP,
                   [&]() { tt.compare_x_on_boundary_2_object()(pt, xcv_bt, xcv_bt_end); });
    t.test_tracing(Tt::COMPARE_X_ON_BOUNDARY_2_OP,
                   [&]() { tt.compare_x_on_boundary_2_object()(xcv_bt, xcv_bt_end, xcv_bt, xcv_bt_end); });
  }
  if constexpr (has_compare_x_near_boundary_2<Tt>::value) {
    t.test_tracing(Tt::COMPARE_X_NEAR_BOUNDARY_2_OP,
                   [&]() { tt.compare_x_near_boundary_2_object()(xcv_bt, xcv_bt, xcv_bt_end); });
  }
  if constexpr (has_compare_y_on_boundary_2<Tt>::value) {
    t.test_tracing(Tt::COMPARE_Y_ON_BOUNDARY_2_OP,
                   [&]() { tt.compare_y_on_boundary_2_object()(pt_lr, pt_lr); });
  }
  if constexpr (has_compare_y_near_boundary_2<Tt>::value) {
    t.test_tracing(Tt::COMPARE_Y_NEAR_BOUNDARY_2_OP,
                   [&]() { tt.compare_y_near_boundary_2_object()(xcv_lr, xcv_lr, xcv_lr_end); });
  }
  if constexpr (has_is_on_x_identification_2<Tt>::value) {
    t.test_tracing(Tt::IS_ON_X_IDENTIFICATION_2_OP, [&]() { tt.is_on_x_identification_2_object()(pt); });
    t.test_tracing(Tt::IS_ON_X_IDENTIFICATION_2_OP, [&]() { tt.is_on_x_identification_2_object()(xcv_seg); });
  }
  if constexpr (has_is_on_y_identification_2<Tt>::value) {
    t.test_tracing(Tt::IS_ON_Y_IDENTIFICATION_2_OP, [&]() { tt.is_on_y_identification_2_object()(pt); });
    t.test_tracing(Tt::IS_ON_Y_IDENTIFICATION_2_OP, [&]() { tt.is_on_y_identification_2_object()(xcv_seg); });
  }
}

void test_tracing_algebraic_traits() {
  using Base_traits = Arr_algebraic_segment_traits_2<CORE::BigInt>;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Tracing_traits, Base_traits>();

  Algebraic_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_bezier_curve_traits() {
  using Base_traits = Bezier_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Bezier_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_circle_segment_traits() {
  using Base_traits = Circle_segment_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Circle_segment_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_conic_traits() {
  using Base_traits = Conic_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Conic_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_consolidated_curve_data_traits() {
  using Base_traits = Consolidated_curve_data_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Consolidated_curve_data_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_curve_data_traits() {
  using Base_traits = Curve_data_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Curve_data_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_geodesic_arc_on_sphere_traits() {
  using Base_traits = Arr_geodesic_arc_on_sphere_traits_2<Exact_predicates_exact_constructions_kernel>;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Tracing_traits, Base_traits>();

  Geodesic_arc_on_sphere_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_linear_traits() {
  using Base_traits = Linear_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Tracing_traits, Base_traits>();

  Linear_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_non_caching_segment_traits() {
  using Base_traits = Non_caching_segment_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Non_caching_segment_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

// void test_tracing_polyline_traits() {
//   using Base_traits = Polyline_traits_decorating_test_types::Base_traits;
//   using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
//   test_type_consistency<Base_traits, Tracing_traits>();

//   Polyline_traits_decorating_test_objects<Tracing_traits> objs;
//   test_tracing_traits(objs);
// }

void test_tracing_rational_funcion_traits() {
  using Base_traits = Rational_function_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Rational_function_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

void test_tracing_segment_traits() {
  using Base_traits = typename Segment_traits_decorating_test_types::Base_traits;
  using Tracing_traits = Arr_tracing_traits_2<Base_traits>;
  test_type_consistency<Base_traits, Tracing_traits>();

  Segment_traits_decorating_test_objects<Tracing_traits> objs;
  test_tracing_traits(objs);
}

int main() {
  test_tracing_algebraic_traits();
  test_tracing_bezier_curve_traits();
  test_tracing_circle_segment_traits();
  test_tracing_conic_traits();
  test_tracing_consolidated_curve_data_traits();
  test_tracing_curve_data_traits();
  test_tracing_geodesic_arc_on_sphere_traits();
  test_tracing_linear_traits();
  test_tracing_non_caching_segment_traits();
  // test_tracing_polyline_traits();
  test_tracing_rational_funcion_traits();
  test_tracing_segment_traits();
  return 0;
}

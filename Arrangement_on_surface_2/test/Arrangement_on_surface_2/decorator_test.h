#ifndef CGAL_DECORATOR_TEST_H
#define CGAL_DECORATOR_TEST_H

#include <CGAL/Arr_has.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE/BigInt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Arr_enums.h>

#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_rational_function_traits_2.h>

using namespace CGAL;
using Epeck = Exact_predicates_exact_constructions_kernel;

/*! Checks the consistency of the presence of functors in two traits-classes.
 */
template <typename T1, typename T2>
void test_type_consistency() {
  using namespace CGAL;

  static_assert(has_compare_x_2<T1>::value == has_compare_x_2<T2>::value, "compare_x_2 presence mismatch");
  static_assert(has_compare_xy_2<T1>::value == has_compare_xy_2<T2>::value, "compare_xy_2 presence mismatch");
  static_assert(has_construct_min_vertex_2<T1>::value == has_construct_min_vertex_2<T2>::value,
                "construct_min_vertex_2 presence mismatch");
  static_assert(has_construct_max_vertex_2<T1>::value == has_construct_max_vertex_2<T2>::value,
                "construct_max_vertex_2 presence mismatch");
  static_assert(has_is_vertical_2<T1>::value == has_is_vertical_2<T2>::value, "is_vertical_2 presence mismatch");
  static_assert(has_compare_y_at_x_2<T1>::value == has_compare_y_at_x_2<T2>::value,
                "compare_y_at_x_2 presence mismatch");
  static_assert(has_equal_2<T1>::value == has_equal_2<T2>::value, "equal_2 presence mismatch");
  static_assert(has_compare_y_at_x_left_2<T1>::value == has_compare_y_at_x_left_2<T2>::value,
                "compare_y_at_x_left_2 presence mismatch");
  static_assert(has_compare_y_at_x_right_2<T1>::value == has_compare_y_at_x_right_2<T2>::value,
                "compare_y_at_x_right_2 presence mismatch");
  static_assert(has_make_x_monotone_2<T1>::value == has_make_x_monotone_2<T2>::value,
                "make_x_monotone_2 presence mismatch");
  static_assert(has_split_2<T1>::value == has_split_2<T2>::value, "split_2 presence mismatch");
  static_assert(has_intersect_2<T1>::value == has_intersect_2<T2>::value, "intersect_2 presence mismatch");
  static_assert(has_are_mergeable_2<T1>::value == has_are_mergeable_2<T2>::value, "are_mergeable_2 presence mismatch");
  static_assert(has_merge_2<T1>::value == has_merge_2<T2>::value, "merge_2 presence mismatch");
  static_assert(has_construct_opposite_2<T1>::value == has_construct_opposite_2<T2>::value,
                "construct_opposite_2 presence mismatch");
  static_assert(has_construct_point_2<T1>::value == has_construct_point_2<T2>::value,
                "construct_point_2 presence mismatch");
  static_assert(has_compare_endpoints_xy_2<T1>::value == has_compare_endpoints_xy_2<T2>::value,
                "compare_endpoints_xy_2 presence mismatch");
  static_assert(has_approximate_2<T1>::value == has_approximate_2<T2>::value, "approximate_2 presence mismatch");
  static_assert(has_approximate_2_point<T1>::value == has_approximate_2_point<T2>::value,
                "approximate_point_2 presence mismatch");
  static_assert(has_approximate_2_xcv<T1>::value == has_approximate_2_xcv<T2>::value,
                "approximate_xcv_2 presence mismatch");
  static_assert(has_approximate_2_xcv_bounds<T1>::value == has_approximate_2_xcv_bounds<T2>::value,
                "approximate_xcv_2_within_bounds presence mismatch");
  static_assert(has_parameter_space_in_x_2<T1>::value == has_parameter_space_in_x_2<T2>::value,
                "parameter_space_in_x_2 presence mismatch");
  static_assert(has_is_on_x_identification_2<T1>::value == has_is_on_x_identification_2<T2>::value,
                "is_on_x_identification_2 presence mismatch");
  static_assert(has_compare_y_on_boundary_2<T1>::value == has_compare_y_on_boundary_2<T2>::value,
                "compare_y_on_boundary_2 presence mismatch");
  static_assert(has_compare_y_near_boundary_2<T1>::value == has_compare_y_near_boundary_2<T2>::value,
                "compare_y_near_boundary_2 presence mismatch");
  static_assert(has_parameter_space_in_y_2<T1>::value == has_parameter_space_in_y_2<T2>::value,
                "parameter_space_in_y_2 presence mismatch");
  static_assert(has_is_on_y_identification_2<T1>::value == has_is_on_y_identification_2<T2>::value,
                "is_on_y_identification_2 presence mismatch");
  static_assert(has_compare_x_on_boundary_2<T1>::value == has_compare_x_on_boundary_2<T2>::value,
                "compare_x_on_boundary_2 presence mismatch");
  static_assert(has_compare_x_near_boundary_2<T1>::value == has_compare_x_near_boundary_2<T2>::value,
                "compare_x_near_boundary_2 presence mismatch");
}

/*! \brief Define a using alias <alias> that is either the existing nested type
 * <namespace>::<identifier> or the provided <fallback_type>.
 */
#define DECL_COND_SNIFAE_TYPE(alias, namespace, identifier, fallback_type)                              \
  template <typename T_, typename = void>                                                               \
  struct conditional_snifae_type_##namespace##_##identifier                                             \
  { using type = fallback_type; };                                                                      \
  template <typename T_>                                                                                \
  struct conditional_snifae_type_##namespace##_##identifier<T_, std::void_t<typename T_::identifier>>   \
  { using type = typename T_::identifier; };                                                            \
  using alias = typename conditional_snifae_type_##namespace##_##identifier<namespace>::type;

template <typename T>
struct Decorating_traits_test_objects {
protected:
  using Point_2 = typename T::Point_2;
  using X_monotone_curve_2 = typename T::X_monotone_curve_2;
  DECL_COND_SNIFAE_TYPE(Curve_2, T, Curve_2, void);

public:
  T traits;
  Point_2 xcv_seg_pt;                                   // a point in the interior of xcv_seg
  X_monotone_curve_2 xcv_seg;                           // x-monotone curve segment
  X_monotone_curve_2 xcv_bt;                            // x-monotone curve approaching bottom/top boundaries (if exists)
  CGAL::Arr_curve_end xcv_bt_end{CGAL::ARR_MIN_END};    // the certain end of xcv_bt that crosses/touches the bottom/top boundary
  X_monotone_curve_2 xcv_lr;                            // x-monotone curve approaching left/right boundaries (if exists)
  CGAL::Arr_curve_end xcv_lr_end{CGAL::ARR_MIN_END};    // the certain end of xcv_lr that crosses/touches the left/right boundary
  Curve_2 cv;                                           // any curve
  Point_2 pt_lr;                                        // a point on the left/right boundary (if exists)
  Point_2 pt_bt;                                        // a point on the bottom/top boundary (if exists)
};

struct Algebraic_traits_decorating_test_types {
  using Base_traits = Arr_algebraic_segment_traits_2<CORE::BigInt>;
};

template <typename T>
struct Algebraic_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Algebraic_traits_decorating_test_objects() {
    using Base_traits = Algebraic_traits_decorating_test_types::Base_traits;
    using Point_2 = typename T::Point_2;
    using Curve_2 = typename T::Curve_2;
    using X_monotone_curve_2 = typename T::X_monotone_curve_2;
    using Polynomial = Base_traits::Polynomial_2;
    using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();

    auto first_x_monotone_curve = [&](const Curve_2& cv) {
      std::vector<Make_x_monotone_result> xcvs;
      base.make_x_monotone_2_object()(cv, std::back_inserter(xcvs));
      return std::get<X_monotone_curve_2>(xcvs.front());
    };

    Polynomial x = CGAL::shift(Polynomial(1), 1, 0);
    Polynomial y = CGAL::shift(Polynomial(1), 1, 1);
    // x^2+y^2-1=0, a unit circle
    this->cv = ctr_cv(CGAL::ipower(x, 4) + CGAL::ipower(y, 2) - 1);
    this->xcv_seg = first_x_monotone_curve(this->cv);
    this->xcv_seg_pt = this->xcv_seg.compare_y_at_x(ctr_pt(0, 1)) == CGAL::EQUAL ? ctr_pt(0, 1) : ctr_pt(0, -1);
    // x^4-y=0
    this->xcv_lr = first_x_monotone_curve(ctr_cv(CGAL::ipower(x, 4) - y));
    // y(x-2)-1=0
    this->xcv_bt = first_x_monotone_curve(ctr_cv((y * (x - 2)) - 1));
    this->xcv_bt_end = ARR_MAX_END;
  }
};

struct Bezier_traits_decorating_test_types {
  using Nt_traits = CORE_algebraic_number_traits;
  using NT = Nt_traits::Rational;
  using Rational = Nt_traits::Rational;
  using Algebraic = Nt_traits::Algebraic;
  using Rat_kernel = Cartesian<Rational>;
  using Alg_kernel = Cartesian<Algebraic>;
  using Rat_point = Rat_kernel::Point_2;
  using Base_traits = Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
};

template <typename T>
struct Bezier_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Bezier_traits_decorating_test_objects() {
    using Point_2 = typename T::Point_2;
    using Curve_2 = typename T::Curve_2;
    using X_monotone_curve_2 = typename T::X_monotone_curve_2;
    using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

    const auto& base = this->traits.traits();

    auto first_x_monotone_curve = [&](const Curve_2& cv) {
      std::vector<Make_x_monotone_result> xcvs;
      base.make_x_monotone_2_object()(cv, std::back_inserter(xcvs));
      return std::get<X_monotone_curve_2>(xcvs.front());
    };

    std::vector<Point_2> pts{{2, 0}, {4, 2}, {0, 4}, {2, 6}};
    this->cv = Curve_2(pts.begin(), pts.end());
    this->xcv_seg = first_x_monotone_curve(this->cv);
    this->xcv_seg_pt = Point_2(2, 3);
  }
};

struct Circle_segment_traits_decorating_test_types {
  using Base_traits = Arr_circle_segment_traits_2<Epeck>;
};

template <typename T>
struct Circle_segment_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Circle_segment_traits_decorating_test_objects() {
    using Base_traits = Circle_segment_traits_decorating_test_types::Base_traits;
    using Rational_point_2 = Base_traits::Rational_point_2;
    using Point_2 = typename T::Point_2;
    using X_monotone_curve_2 = typename T::X_monotone_curve_2;
    using Circle_2 = Base_traits::Rational_circle_2;

    auto circ = Circle_2(Rational_point_2(0, 0), 4);
    this->cv = circ;
    this->xcv_seg = X_monotone_curve_2(circ, Point_2(-2, 0), Point_2(2, 0), CLOCKWISE);
    this->xcv_seg_pt = Point_2(0, 2);
  }
};

struct Conic_traits_decorating_test_types {
  using Nt_traits = CORE_algebraic_number_traits;
  using Rational = Nt_traits::Rational;
  using Algebraic = Nt_traits::Algebraic;
  using Rat_kernel = Cartesian<Rational>;
  using Alg_kernel = Cartesian<Algebraic>;
  using Base_traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
};

template <typename T>
struct Conic_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Conic_traits_decorating_test_objects() {
    using Point_2 = typename T::Point_2;
    using Rational = typename Conic_traits_decorating_test_types::Rational;

    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();

    this->cv = ctr_cv(0, 0, 1, 0, 0, -1, COUNTERCLOCKWISE, Point_2(Rational(1, 4), 4), Point_2(2, Rational(1, 2)));
    this->xcv_seg = ctr_xcv(this->cv);
    this->xcv_seg_pt = Point_2(Rational(1, 1), 1);
  }
};

struct Consolidated_curve_data_traits_decorating_test_types {
  using Segment_traits = Arr_segment_traits_2<Epeck>;
  using Base_traits = Arr_consolidated_curve_data_traits_2<Segment_traits, int>;
};

template <typename T>
struct Consolidated_curve_data_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Consolidated_curve_data_traits_decorating_test_objects() {
    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();
    this->xcv_seg = ctr_xcv(ctr_pt(0, 0), ctr_pt(1, 1), 0);
    this->xcv_seg_pt = ctr_pt(0.5, 0.5);
    this->cv = ctr_cv(ctr_pt(0, 0), ctr_pt(1, 1));
  }
};

struct Curve_data_traits_decorating_test_types {
  using Segment_traits = Arr_segment_traits_2<Epeck>;
  using Base_traits = Arr_curve_data_traits_2<Segment_traits, int>;
};

template <typename T>
struct Curve_data_traits_decorating_test_objects :
    public Decorating_traits_test_objects<T> {
  Curve_data_traits_decorating_test_objects() {
    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();
    this->xcv_seg = ctr_xcv(ctr_pt(0, 0), ctr_pt(1, 1), 0);
    this->xcv_seg_pt = ctr_pt(0.5, 0.5);
    this->cv = ctr_cv(ctr_pt(0, 0), ctr_pt(1, 1));
  }
};

struct Geodesic_arc_on_sphere_traits_decorating_test_types {
  using Base_traits = Arr_geodesic_arc_on_sphere_traits_2<Epeck>;
};

template <typename T>
struct Geodesic_arc_on_sphere_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Geodesic_arc_on_sphere_traits_decorating_test_objects() {
    const auto& base = this->traits.traits();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();
    this->xcv_seg = ctr_xcv(ctr_pt(1, 0, 0), ctr_pt(0, 1, 0));
    this->xcv_seg_pt = ctr_pt(0.5, 0.5, 0);
    this->cv = ctr_cv(ctr_pt(1, 0, 0), ctr_pt(0, 0, 1));
    this->xcv_bt = ctr_xcv(ctr_pt(0, 1, 0), ctr_pt(0, 0, 1));
    this->xcv_bt_end = ARR_MAX_END;
    this->xcv_lr = ctr_xcv(ctr_pt(-1, 0, -1), ctr_pt(-1, 0, 1));
    this->pt_lr = ctr_pt(-1, 0, 0);
    this->pt_bt = ctr_pt(0, 0, 1);
  }
};

struct Linear_traits_decorating_test_types {
  using Base_traits = Arr_linear_traits_2<Epeck>;
};

template <typename T>
struct Linear_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Linear_traits_decorating_test_objects() {
    using Base_traits = Linear_traits_decorating_test_types::Base_traits;
    using Point_2 = typename T::Point_2;
    using X_monotone_curve_2 = typename T::X_monotone_curve_2;

    this->xcv_seg_pt = Point_2(0.5, 0.5);
    this->xcv_seg = X_monotone_curve_2(Point_2(0, 0), Point_2(1, 1));
    this->xcv_lr = Base_traits::Line_2(Point_2(0, 0), Point_2(1, 0));
    this->xcv_bt = Base_traits::Line_2(Point_2(0, 0), Point_2(0, 1));
    this->cv = this->xcv_seg;
  }
};

struct Non_caching_segment_traits_decorating_test_types {
  using Base_traits = Arr_non_caching_segment_traits_2<Epeck>;
};

template <typename T>
struct Non_caching_segment_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Non_caching_segment_traits_decorating_test_objects() {
    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();

    this->xcv_seg = ctr_xcv(ctr_pt(0, 0), ctr_pt(1, 1));
    this->xcv_seg_pt = ctr_pt(0.5, 0.5);
    this->cv = ctr_cv(ctr_pt(0, 0), ctr_pt(1, 1));
  }
};

struct Polyline_traits_decorating_test_types {
  using Base_traits = Arr_polyline_traits_2<>;
};

template <typename T>
struct Polyline_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Polyline_traits_decorating_test_objects() {
    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();

    this->xcv_seg = ctr_xcv(ctr_pt(0, 0), ctr_pt(1, 1));
    this->xcv_seg_pt = ctr_pt(0.5, 0.5);
    this->cv = ctr_cv(ctr_pt(0, 0), ctr_pt(1, 1));
  }
};

struct Rational_function_traits_decorating_test_types {
  using Number_type = CORE::BigInt;
  using AK1 = CGAL::Algebraic_kernel_d_1<Number_type>;
  using Base_traits = CGAL::Arr_rational_function_traits_2<AK1>;
};

template <typename T>
struct Rational_function_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Rational_function_traits_decorating_test_objects() {
    using Base_traits = typename Rational_function_traits_decorating_test_types::Base_traits;
    using AK1 = typename Rational_function_traits_decorating_test_types::AK1;
    using Alg_real = typename AK1::Algebraic_real_1;
    using Polynomial = typename AK1::Polynomial_1;
    using Rational = typename Base_traits::Rational;
    using Bound = typename AK1::Bound;

    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();

    Polynomial x = CGAL::shift(Polynomial(1), 1);
    Polynomial P1 = CGAL::ipower(x, 4) - 6 * x * x + 8;
    this->xcv_seg = ctr_xcv(P1, Polynomial(1), Alg_real(Bound(-2)), Alg_real(Bound(2)));
    this->xcv_seg_pt = ctr_pt(Rational(0), P1.evaluate(0));
    this->cv = ctr_cv(P1, Polynomial(1));
    this->xcv_lr = ctr_xcv(P1, Polynomial(1));
    this->xcv_lr_end = ARR_MIN_END;
    this->xcv_bt = ctr_xcv(Polynomial(1), x, Alg_real(Bound(0)), Alg_real(Bound(2)));
    this->xcv_bt_end = ARR_MIN_END;
  }
};

struct Segment_traits_decorating_test_types {
  using Base_traits = Arr_segment_traits_2<Epeck>;
};

template <typename T>
struct Segment_traits_decorating_test_objects : public Decorating_traits_test_objects<T> {
  Segment_traits_decorating_test_objects() {
    const auto& base = this->traits.traits();
    auto ctr_cv = base.construct_curve_2_object();
    auto ctr_pt = base.construct_point_2_object();
    auto ctr_xcv = base.construct_x_monotone_curve_2_object();

    this->xcv_seg = ctr_xcv(ctr_pt(0, 0), ctr_pt(1, 1));
    this->xcv_seg_pt = ctr_pt(0.5, 0.5);
    this->cv = ctr_cv(ctr_pt(0, 0), ctr_pt(1, 1));
  }
};

template <typename, typename = std::void_t<>>
struct is_basic_traits : std::false_type {};

template <typename T>
struct is_basic_traits<T, std::void_t<
  typename T::Point_2,
  typename T::X_monotone_curve_2,
  typename T::Has_left_category,
  typename T::Left_side_category,
  typename T::Bottom_side_category,
  typename T::Top_side_category,
  typename T::Right_side_category,
  typename T::Compare_x_2,
  typename T::Compare_xy_2,
  typename T::Construct_min_vertex_2,
  typename T::Construct_max_vertex_2,
  typename T::Is_vertical_2,
  typename T::Compare_y_at_x_2,
  typename T::Compare_y_at_x_left_2,
  typename T::Compare_y_at_x_right_2,
  typename T::Equal_2
>> : std::true_type {};

// Detect whether T models ArrangementBasicTraits_2
template <typename T>
constexpr bool is_basic_traits_v = is_basic_traits<T>::value;

template <typename T>
struct unwrap { using type = T; };

template <template <typename> class Wrapped, typename T>
struct unwrap<Wrapped<T>>
{ using type = std::conditional_t<is_basic_traits_v<T>, typename unwrap<T>::type, Wrapped<T>>; };

// Recursively extract the underlying basic traits type
template <typename T>
using unwrap_t = typename unwrap<T>::type;

// Detect whether Gt is or derives from Arr_conic_traits_2<*, *, *>
template <typename Gt>
struct is_or_derived_from_conic_traits {
private:
  template <typename RatKernel, typename AlgKernel, typename NtTraits>
  static std::true_type test(const Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>*);

  static std::false_type test(...);

public:
  static constexpr bool value = decltype(test(static_cast<const Gt*>(nullptr)))::value;
};

template <typename Gt>
inline constexpr bool is_or_derived_from_conic_traits_v = is_or_derived_from_conic_traits<Gt>::value;


#endif

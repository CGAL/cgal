#ifndef CGAL_DRAW_AOS_TYPE_UTILS_H
#define CGAL_DRAW_AOS_TYPE_UTILS_H
#include "CGAL/Arr_polycurve_traits_2.h"
#include "CGAL/Arr_polyline_traits_2.h"
#include "CGAL/Cartesian.h"
#include <type_traits>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>

namespace CGAL {
namespace draw_aos {

enum class Side_of_boundary {
  Top = 0,
  Left = 1,
  Bottom = 2,
  Right = 3,
  None = -1,
};

template <typename, typename = std::void_t<>>
struct has_construct_x_monotone_curve_2 : std::false_type
{};

template <typename T>
struct has_construct_x_monotone_curve_2<T, std::void_t<typename T::Construct_x_monotone_curve_2>> : std::true_type
{};

template <typename, typename = std::void_t<>>
struct has_approximate_2_object : std::false_type
{};

// Specialization: detection succeeds if decltype(T::approximate_2_object()) is valid
template <typename T>
struct has_approximate_2_object<T, std::void_t<decltype(std::declval<T>().approximate_2_object())>> : std::true_type
{};

// Convenience variable
template <typename T>
inline constexpr bool has_approximate_2_object_v = has_approximate_2_object<T>::value;

// Primary templates: detection fails by default
// Does a class have operator()(const Point&)?
template <typename, typename, typename = std::void_t<>>
struct has_operator_point : std::false_type
{};

// Specialization: detection succeeds if decltype works out
template <typename T, typename A>
struct has_operator_point<T, A, std::void_t<decltype(std::declval<A>()(std::declval<const typename T::Point_2&>()))>>
    : std::true_type
{};

// Convenience variable
template <typename T, typename A>
inline constexpr bool has_operator_point_v = has_operator_point<T, A>::value;

// Primary templates: detection fails by default
// Does a class have operator()(const X_monotone_curve&)?
template <typename, typename, typename, typename = std::void_t<>>
struct has_operator_xcv : std::false_type
{};

template <typename T, typename A, typename O>
struct has_operator_xcv<T,
                        A,
                        O,
                        std::void_t<decltype(std::declval<A&>()(std::declval<const typename T::X_monotone_curve_2&>(),
                                                                std::declval<double>(),
                                                                std::declval<O>(),
                                                                std::declval<bool>()))>> : std::true_type
{};

// Convenience variable
template <typename T, typename A>
constexpr bool has_operator_xcv_v = has_operator_xcv<T, A, void*>::value;

template <typename GeomTraits>
struct Traits_adaptor_base
{
public:
  using Geom_traits = GeomTraits;
  using Point_2 = typename Geom_traits::Point_2;
  using X_monotone_curve_2 = typename Geom_traits::X_monotone_curve_2;
  using Intersect_2 = typename Geom_traits::Intersect_2;
  using Is_vertical_2 = typename Geom_traits::Is_vertical_2;
  using Compare_xy_2 = typename Geom_traits::Compare_xy_2;
};

template <typename GeomTraits>
struct Traits_adaptor;

template <typename Kernel>
struct Traits_adaptor<Arr_segment_traits_2<Kernel>> : public Traits_adaptor_base<Arr_segment_traits_2<Kernel>>
{
private:
  using Geom_traits = Arr_segment_traits_2<Kernel>;

public:
  constexpr static bool Has_unbounded_curves = false;
  constexpr static double Approximation_sizing_factor = 1.0;
  using FT = typename Kernel::FT;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using Approximate_number_type = typename Geom_traits::Approximate_number_type;
  using Approximate_kernel = typename Geom_traits::Approximate_kernel;
  using Approximate_point_2 = typename Geom_traits::Approximate_point_2;
};

template <typename SegmentTraits>
struct Traits_adaptor<Arr_polyline_traits_2<SegmentTraits>>
    : public Traits_adaptor_base<Arr_polyline_traits_2<SegmentTraits>>
{
private:
  using Geom_traits = Arr_polyline_traits_2<SegmentTraits>;
  using Sub_traits = SegmentTraits;
  using Adapted_sub_traits = Traits_adaptor<Sub_traits>;

public:
  constexpr static bool Has_unbounded_curves = false;
  constexpr static double Approximation_sizing_factor = 1.0;
  using FT = typename Adapted_sub_traits::FT;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using Approximate_number_type = typename Adapted_sub_traits::Approximate_number_type;
  using Approximate_kernel = typename Adapted_sub_traits::Approximate_kernel;
  using Approximate_point_2 = typename Adapted_sub_traits::Approximate_point_2;
};

template <typename SubcurveTraits>
struct Traits_adaptor<Arr_polycurve_traits_2<SubcurveTraits>>
    : public Traits_adaptor_base<Arr_polycurve_traits_2<SubcurveTraits>>
{
private:
  using Sub_traits = SubcurveTraits;
  using Geom_traits = Arr_polycurve_traits_2<Sub_traits>;
  using Adapted_sub_traits = Traits_adaptor<Sub_traits>;

public:
  constexpr static bool Has_unbounded_curves = false;
  constexpr static double Approximation_sizing_factor = 1.0;
  using FT = typename Adapted_sub_traits::FT;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using Approximate_number_type = typename Adapted_sub_traits::Approximate_number_type;
  using Approximate_kernel = typename Adapted_sub_traits::Approximate_kernel;
  using Approximate_point_2 = typename Adapted_sub_traits::Approximate_point_2;
};

template <typename Kernel>
struct Traits_adaptor<Arr_linear_traits_2<Kernel>> : public Traits_adaptor_base<Arr_linear_traits_2<Kernel>>
{
private:
  using Geom_traits = Arr_segment_traits_2<Kernel>;

public:
  constexpr static bool Has_unbounded_curves = true;
  constexpr static double Approximation_sizing_factor = 0.0;
  using FT = typename Kernel::FT;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using Approximate_number_type = typename Geom_traits::Approximate_number_type;
  using Approximate_kernel = typename Geom_traits::Approximate_kernel;
  using Approximate_point_2 = typename Geom_traits::Approximate_point_2;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
struct Traits_adaptor<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>
    : public Traits_adaptor_base<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>
{
private:
  using Geom_traits = Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>;

public:
  constexpr static bool Has_unbounded_curves = false;
  constexpr static double Approximation_sizing_factor = 1.0;
  using FT = typename AlgKernel::FT;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using Approximate_number_type = typename Geom_traits::Approximate_number_type;
  using Approximate_kernel = typename Geom_traits::Approximate_kernel;
  using Approximate_point_2 = typename Geom_traits::Approximate_point_2;
};

template <typename Kernel>
struct Traits_adaptor<Arr_circle_segment_traits_2<Kernel>>
    : public Traits_adaptor_base<Arr_circle_segment_traits_2<Kernel>>
{
private:
  using Geom_traits = Arr_circle_segment_traits_2<Kernel>;
  using Base = Traits_adaptor_base<Geom_traits>;

public:
  constexpr static bool Has_unbounded_curves = false;
  constexpr static double Approximation_sizing_factor = 0.5;
  using FT = typename Base::Point_2::CoordNT;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using Approximate_number_type = typename Geom_traits::Approximate_number_type;
  using Approximate_kernel = typename Geom_traits::Approximate_kernel;
  using Approximate_point_2 = typename Geom_traits::Approximate_point_2;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
struct Traits_adaptor<Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>
    : public Traits_adaptor_base<Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>
{
  static_assert(false, "Approximate_2 not yet modeled by this geometry traits class.");
};

template <typename Kernel>
struct Traits_adaptor<Arr_circular_line_arc_traits_2<Kernel>>
    : public Traits_adaptor_base<Arr_circular_line_arc_traits_2<Kernel>>
{
  static_assert(false, "Approximate_2 not yet modeled by this geometry traits class.");
};

template <typename Kernel>
struct Traits_adaptor<Arr_rational_function_traits_2<Kernel>>
    : public Traits_adaptor_base<Arr_rational_function_traits_2<Kernel>>
{
private:
  using Geom_traits = Arr_rational_function_traits_2<Kernel>;

public:
  constexpr static bool Has_unbounded_curves = true;
  constexpr static double Approximation_sizing_factor = 5.0;
  using FT = typename Geom_traits::Algebraic_real_1;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  // Currently, Approximate_number_type is defined as Bound in Arr_rational_function_traits_2,
  // And there's no Approximate_kernel defined.
  using Approximate_number_type = double;
  using Approximate_kernel = Cartesian<double>;
  using Approximate_point_2 = typename Approximate_kernel::Point_2;
};

template <typename GeomTraits>
class Construct_coordinate
{
  using FT = typename Traits_adaptor<GeomTraits>::FT;

public:
  FT operator()(double val) const { return FT(val); }
};

template <typename Kernel>
class Construct_coordinate<Arr_rational_function_traits_2<Kernel>>
{
  using FT = typename Traits_adaptor<Arr_rational_function_traits_2<Kernel>>::FT;
  using Bound = typename Kernel::Bound;

public:
  FT operator()(double val) const { return FT(Bound(val)); }
};

template <typename GeomTraits>
class Arr_approximation_geometry_traits
{
  using Adapted_traits = Traits_adaptor<GeomTraits>;

public:
  using Approx_point = typename Adapted_traits::Approximate_point_2;
  using Approx_nt = typename Adapted_traits::Approximate_number_type;
  using Approx_kernel = typename Adapted_traits::Approximate_kernel;
  using Point_geom = Approx_point;
  using Apporx_point_vec = std::vector<Point_geom>;
  using Polyline_geom = Apporx_point_vec;
  using Triangle = std::array<std::size_t, 3>;
  using Triangle_vec = std::vector<Triangle>;
  struct Triangulated_face
  {
    Apporx_point_vec points;
    Triangle_vec triangles;
  };
};

} // namespace draw_aos
} // namespace CGAL

#endif // CGAL_DRAW_AOS_TYPE_UTILS_H
#ifndef CGAL_DRAW_AOS_TYPE_UTILS_H
#define CGAL_DRAW_AOS_TYPE_UTILS_H
#include <vector>
#include <type_traits>

#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

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

/*!
 */
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

template <typename T>
constexpr bool has_approximate_traits_v =
    has_approximate_2_object_v<T> && has_operator_point_v<T, typename T::Approximate_2> &&
    has_operator_xcv_v<T, typename T::Approximate_2>;

// Detect whether T is or derives from Arr_geodesic_arc_on_sphere_traits_2<*, *, *>
template <typename T>
struct is_or_derived_from_agas
{
private:
  template <typename Kernel_, int AtanX, int AtanY>
  static std::true_type test(const Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY>*);

  static std::false_type test(...);

public:
  static constexpr bool value = decltype(test(static_cast<const T*>(nullptr)))::value;
};

template <typename T>
inline constexpr bool is_or_derived_from_agas_v = is_or_derived_from_agas<T>::value;

// Detect whether T is or derives from a geometry traits on curved surfaces
template <typename T>
inline constexpr bool is_or_derived_from_curved_surf_traits = is_or_derived_from_agas_v<T>;

/*!
 */
template <typename GeomTraits>
class Arr_approximate_traits
{
  using Geom_traits = GeomTraits;

public:
  using Approx_point = typename Geom_traits::Approximate_point_2;
  using Approx_nt = typename Geom_traits::Approximate_number_type;
  using Approx_kernel = typename Geom_traits::Approximate_kernel;
  using Approx_proj_point = typename Approx_kernel::Point_2;
  using Point_geom = Approx_proj_point;
  using Point_geom_vec = std::vector<Point_geom>;
  using Polyline_geom = Point_geom_vec;
  using Triangle = std::array<std::size_t, 3>;
  using Triangle_vec = std::vector<Triangle>;
  struct Triangulated_face
  {
    Point_geom_vec points;
    Triangle_vec triangles;
  };
};

} // namespace draw_aos
} // namespace CGAL

#endif // CGAL_DRAW_AOS_TYPE_UTILS_H

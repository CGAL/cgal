#ifndef CGAL_DRAW_AOS_ARR_PROJECTION_H
#define CGAL_DRAW_AOS_ARR_PROJECTION_H
#include <cmath>

#include <CGAL/number_type_config.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/**
 * @brief class handling projection between 2D parameter space and coordinate space
 *
 * @tparam GeomTraits
 */
template <typename GeomTraits>
class Arr_projector
{
  using Geom_traits = GeomTraits;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approx_proj_point = typename Approx_traits::Approx_proj_point;

public:
  Arr_projector(const GeomTraits& traits)
      : m_traits(traits) {}

  Approx_proj_point project(Approx_point pt) const { return pt; }

  Approx_point unproject(Approx_proj_point pt) const { return pt; }

private:
  const GeomTraits& m_traits;
};

/**
 * @brief Projector specialization for geodesic arc on sphere traits.
 *
 * The projection process is essentially a conversion between spherical coordinates and right-handed Cartesian
 * coordinates. Sphercial coordinates are represented as azimuth and polar angle. Note that the polar angle starts from
 * north pole.
 *
 * @tparam Kernel
 * @tparam atanX
 * @tparam atanY
 */
template <typename Kernel, int atanX, int atanY>
class Arr_projector<Arr_geodesic_arc_on_sphere_traits_2<Kernel, atanX, atanY>>
{
  using Geom_traits = Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Approx_proj_point = typename Approx_traits::Approx_proj_point;

public:
  Arr_projector(const Geom_traits& traits)
      : m_traits(traits) {}

  Approx_proj_point project(Approx_point point) const {
    if(point.location() == Approx_point::MAX_BOUNDARY_LOC) return Approx_proj_point(0, CGAL_PI);
    if(point.location() == Approx_point::MIN_BOUNDARY_LOC) return Approx_proj_point(0, 0);
    Approx_nt azimuth_from_id =
        std::fmod(std::atan2(point.dy(), point.dx()) - std::atan2(atanY, atanX) + 2 * CGAL_PI, 2 * CGAL_PI);
    return Approx_proj_point(azimuth_from_id, std::acos(-point.dz()));
  }

  Approx_point unproject(Approx_proj_point point) const {
    using Direction_3 = typename Geom_traits::Approximate_kernel::Direction_3;

    Approx_nt polar = point.y();
    if(point.y() == CGAL_PI) return Approx_point(Direction_3(0, 0, 1), Approx_point::MAX_BOUNDARY_LOC);
    if(point.y() == 0) return Approx_point(Direction_3(0, 0, -1), Approx_point::MIN_BOUNDARY_LOC);
    Approx_nt azimuth = point.x() + std::atan2(atanY, atanX);
    return Approx_point(
        Direction_3(std::sin(polar) * std::cos(azimuth), std::sin(polar) * std::sin(azimuth), -std::cos(polar)),
        azimuth == 0 ? Approx_point::MID_BOUNDARY_LOC : Approx_point::NO_BOUNDARY_LOC);
  }

private:
  const Geom_traits& m_traits;
};

} // namespace draw_aos
} // namespace CGAL

#endif
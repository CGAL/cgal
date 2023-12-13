// Copyright (c) 2023 Tel-Aviv University (Israel).
namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Ref
 *
 * \anchor arr_refarr_obs
 *
 * `Arr_observer<Arrangement_2>` is an alias for
 * Aos_observer<Arrangement_on_surface_2>`,
 * where `Arrangement_2` derives from `Arrangement_on_surface_2` and the latter
 * is an instance of the template
 * `CGAL::Arrangement_on_surface_2<GeometryTraits, TopologyTraits>`.
 */

template <typename Arrangement_>
using Arr_observer = typename Arrangement_::Observer;

} // namespace CGAL

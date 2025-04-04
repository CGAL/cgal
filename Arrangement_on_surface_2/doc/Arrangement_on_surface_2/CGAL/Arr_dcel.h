namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2DCEL
 *
 * The \dcel class used by the `Arrangement_2`,
 * `Arr_bounded_planar_topology_traits_2`, `Arr_unb_planar_topology_traits_2`
 * class templates and other templates.  It is parameterized by a geometry
 * traits type and optionally by a vertex, halfedge, or face types. By default,
 * the `Arr_dcel` class template uses the \link
 * AosBasicTraits_2::Point_2 `Point_2`\endlink and \link
 * AosBasicTraits_2::X_monotone_curve_2 `X_monotone_curve_2`\endlink
 * types nested in the traits type to instantiate the vertex and base halfedge
 * types, respectively. Thus, by default the \dcel only stores the topological
 * incidence relations and the geometric data attached to vertices and
 * edges. Any one of the vertex, halfedge, or face types can be
 * overridden. Notice that if the vertex and halfedge types are overridden, the
 * traits type is ignored.
 *
 * \cgalModels{AosDcelWithRebind}
 *
 * \tparam Traits a geometry traits type, which is a model of the
 *                `AosBasicTraits_2` concept.
 * \tparam V the vertex type, which is a model of the `AosDcelVertex`
 *           concept.
 * \tparam H the halfedge type, which is a model of the
 *            `AosDcelHalfedge` concept.
 * \tparam F the face type, which is a model of the `AosDcelFace`
 *           concept.
 *
 * \sa `Arr_dcel_base<V, H, F>`
 */
template <typename Traits,
          typename V = Arr_vertex_base<typename Traits::Point_2>,
          typename H = Arr_halfedge_base<typename Traits::X_monotone_curve_2>,
          typename F = Arr_face_base>
class Arr_dcel : public Arr_dcel_base<V, H, F> {};

} /* end namespace CGAL */

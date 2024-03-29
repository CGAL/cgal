namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2DCEL
 *
 * The default \dcel class used by the `Arrangement_2`,
 * `Arr_bounded_planar_topology_traits_2`, `Arr_unb_planar_topology_traits_2`
 * class templates and other templates.  It is parameterized by a geometry
 * traits type. It uses the \link ArrangementBasicTraits_2::Point_2
 * `Point_2`\endlink and \link ArrangementBasicTraits_2::X_monotone_curve_2
 * `X_monotone_curve_2`\endlink types nested in the traits type to instantiate
 * the vertex and base halfedge types, respectively. Thus, by default the \dcel
 * only stores the topological incidence relations and the geometric data
 * attached to vertices and edges.
 *
 * \cgalModels{ArrangementDcelWithRebind}
 *
 * \tparam Traits a geometry traits type, which is a model of the
 *                `ArrangementBasicTraits_2` concept.
 *
 * \sa `Arr_dcel<Traits, V, H, F>`
 * \sa `Arr_dcel_base<V, H, F>`
 */
template <typename Traits> using Arr_default_dcel = Arr_dcel<Traits>;

} /* end namespace CGAL */

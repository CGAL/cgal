/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

The concept `EdgeIntersectionOracle_3` describes the requirements for the edge-isosurface
intersection oracle template parameter of the domain classes
`CGAL::Isosurfacing::Marching_cubes_domain_3` and
`CGAL::Isosurfacing::Dual_contouring_domain_3`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Dichotomy_edge_intersection}
\cgalHasModels{CGAL::Isosurfacing::Linear_interpolation_edge_intersection}
\cgalHasModelsEnd
*/
class EdgeIntersectionOracle_3
{
public:
  /*!
   * \brief computes the intersection point between an edge and the isosurface.
   *
   * \tparam Domain must be a model of `IsosurfacingDomain_3`
   *
   * \param p_0 the geometric position of the first vertex of the edge
   * \param p_1 the geometric position of the second vertex of the edge
   * \param val_0 the value at the first vertex of the edge
   * \param val_1 the value at the second vertex of the edge
   * \param domain the isosurfacing domain
   * \param isovalue the isovalue defining the isosurfacing with which we seek an intersection
   * \param p the intersection point, if it exists
   *
   * \return `true` if the intersection point exists, `false` otherwise
   */
  template <typename Domain>
  bool operator()(const typename Domain::Geom_traits::Point_3& p_0,
                  const typename Domain::Geom_traits::Point_3& p_1,
                  const typename Domain::Geom_traits::FT val_0,
                  const typename Domain::Geom_traits::FT val_1,
                  const Domain& domain,
                  const typename Domain::Geom_traits::FT isovalue,
                  typename Domain::Geom_traits::Point_3& p) const;
};

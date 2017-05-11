/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * Models of the concept `ArrangementSphericalBoundaryTraits_2` handle curves on
 * a sphere or a surface that is topological equivalent to a sphere. The sphere
 * is oriented in such a way that the boundary of the rectangular parameter
 * space, the sphere is the mapping of which, is identified on the left and
 * right sides and contracted at the top and bottom sides. In other words,
 *
 * \cgalRefines `ArrangementBasicTraits_2`
 * \cgalRefines `ArrangementIdentifiedLeftTraits_2`
 * \cgalRefines `ArrangementIdentifiedRightTraits_2`
 * \cgalRefines `ArrangementContractedBottomTraits_2`
 * \cgalRefines `ArrangementContractedTopTraits_2`
 * \cgalRefines
 *
 * \cgalHasModel `CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>`
 *
 * \sa `ArrangementOpenBoundaryTraits_2`
 * \sa `ArrangementBasicTraits_2`
 * \sa `ArrangementIdentifiedLeftTraits_2`
 * \sa `ArrangementIdentifiedRightTraits_2`
 * \sa `ArrangementContractedBottomTraits_2`
 * \sa `ArrangementContractedTopTraits_2`
 * \sa `ArrangementHorizontalSideTraits_2`
 * \sa `ArrangementVerticalSideTraits_2`
 */

class ArrangementSphericalBoundaryTraits_2 {
public:

  /// \name Categories
  /// @{
  /// @}

  /// \name Functor Types
  /// @{
  /// @}

}; /* end ArrangementSphericalBoundaryTraits_2 */

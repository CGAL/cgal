/*!
 * \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementContractedBottomTraits_2` must be used
 * when the parameter space of the surface, the arrangement is embedded on, is
 * contracted on the bottom side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the bottom boundary side when it is contracted.
 *
 * \cgalRefines{ArrangementBottomSideTraits_2}
 *
 * \sa `ArrangementContractedLeftTraits_2`,
 *     `ArrangementContractedRightTraits_2`,
 *     `ArrangementContractedTopTraits_2`,
 *     `ArrangementClosedBottomTraits_2`,
 *     `ArrangementContractedBottomTraits_2`, and
 *     `ArrangementIdentifiedHorizontalTraits_2`
 */

class ArrangementContractedBottomTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_contracted_side_tag`.
  typedef unspecified_type Bottom_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// \name Accessing Functor Objects
  /// @{
  /// @}
}

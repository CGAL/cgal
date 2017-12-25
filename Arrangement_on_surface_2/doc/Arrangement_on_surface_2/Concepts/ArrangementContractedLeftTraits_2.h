/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementContractedLeftTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * contracted on the left side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the left boundary side when it is contracted.
 *
 * \cgalRefines `ArrangementLeftSideTraits_2`
 *
 * \sa `ArrangementContractedRightTraits_2`,
 *     `ArrangementContractedBottomTraits_2`,
 *     `ArrangementContractedTopTraits_2`,
 *     `ArrangementOpenLeftTraits_2`,
 *     `ArrangementClosedLeftTraits_2`, and
 *     `ArrangementIdentifiedLeftTraits_2`
 */

class ArrangementContractedLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_contracted_side_tag`.
  typedef unspecified_type Left_side_category;
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

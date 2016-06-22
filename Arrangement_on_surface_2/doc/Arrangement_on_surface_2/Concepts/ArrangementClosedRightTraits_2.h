/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementClosedRightTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is closed on
 * the right side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * right boundary side when it is closed.
 *
 * \cgalRefines `ArrangementRightSideTraits_2`
 *
 * \sa `ArrangementClosedLeftTraits_2`,
 *     `ArrangementClosedBottomTraits_2`,
 *     `ArrangementClosedTopTraits_2`,
 *     `ArrangementOpenRightTraits_2`,
 *     `ArrangementContractedRightTraits_2`, and
 *     `ArrangementIdentifiedRightTraits_2`,
 */

class ArrangementClosedRightTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_closed_side_tag`.
  typedef unspecified_type Right_side_category;
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

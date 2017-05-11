/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementClosedTopTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is closed on
 * the top side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * top boundary side when it is closed.
 *
 * \cgalRefines `ArrangementTopSideTraits_2`
 *
 * \sa `ArrangementClosedLeftTraits_2`,
 *     `ArrangementClosedRightTraits_2`,
 *     `ArrangementClosedBottomTraits_2`,
 *     `ArrangementOpenTopTraits_2`,
 *     `ArrangementContractedTopTraits_2`, and
 *     `ArrangementIdentifiedTopTraits_2`
 */

class ArrangementClosedTopTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_closed_side_tag`.
  typedef unspecified_type Top_side_category;
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

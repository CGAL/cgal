/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementOpenLeftTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the left side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * left boundary side when it is open.
 *
 * \cgalRefines `ArrangementLeftSideTraits_2`
 *
 * \sa `ArrangementOpenRightTraits_2`,
 *     `ArrangementOpenBottomTraits_2`,
 *     `ArrangementOpenTopTraits_2`,
 *     `ArrangementClosedLeftTraits_2`,
 *     `ArrangementContractedLeftTraits_2`, and
 *     `ArrangementIdentifiedLeftTraits_2`
 */

class ArrangementOpenLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_open_side_tag`.
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

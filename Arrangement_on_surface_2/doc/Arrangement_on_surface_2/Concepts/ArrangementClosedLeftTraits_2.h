/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementClosedLeftTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is closed on
 * the left side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * left boundary side when it is closed.

 * \cgalRefines `ArrangementLeftSideTraits_2`

\sa `ArrangementClosedRightTraits_2`,
    `ArrangementClosedBottomTraits_2`,
    `ArrangementClosedTopTraits_2`,
    `ArrangementOpenLeftTraits_2`,
    `ArrangementContractedLeftTraits_2`, and
    `ArrangementIdentifiedLeftTraits_2`,
*/

class ArrangementClosedLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  /*! Must be convertible to `Arr_closed_side_tag`.
  */
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

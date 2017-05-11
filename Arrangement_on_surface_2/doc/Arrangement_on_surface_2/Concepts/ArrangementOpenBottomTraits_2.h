/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementOpenBottomTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the bottom side and curves inserted into the arrangement are expected to
 * reach this boundary side. A model of this concept can handle curves that
 * reach the bottom boundary side when it is open.
 *
 * \cgalRefines `ArrangementBottomSideTraits_2`
 *
 * \sa `ArrangementOpenLeftTraits_2`,
 *     `ArrangementOpenRightTraits_2`,
 *     `ArrangementOpenTopTraits_2`,
 *     `ArrangementClosedBottomTraits_2`,
 *     `ArrangementContractedBottomTraits_2`, and
 *     `ArrangementIdentifiedBottomTraits_2`
 */

class ArrangementOpenBottomTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_open_side_tag`.
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

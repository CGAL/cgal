/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementIdentifiedBottomTraits_2` must be used
 * when the parameter space of the surface, the arrangement is embedded on, is
 * identified on the bottom side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the bottom boundary side when it is identified.
 *
 * \cgalRefines `ArrangementBottomSideTraits_2`
 *
 * \sa `ArrangementIdentifiedRightTraits_2`,
 *     `ArrangementIdentifiedBottomTraits_2`,
 *     `ArrangementIdentifiedTopTraits_2`,
 *     `ArrangementOpenBottomTraits_2`,
 *     `ArrangementClosedBottomTraits_2`, and
 *     `ArrangementConstractedBottomTraits_2`
 */

class ArrangementIdentifiedBottomTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_identified_side_tag`.
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

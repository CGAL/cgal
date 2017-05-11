/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementIdentifiedRightTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * identified on the right side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the right boundary side when it is identified.
 *
 * \cgalRefines `ArrangementRightSideTraits_2`
 *
 * \sa `ArrangementIdentifiedLeftTraits_2`,
 *     `ArrangementIdentifiedBottomTraits_2`,
 *     `ArrangementIdentifiedTopTraits_2`,
 *     `ArrangementOpenRightTraits_2`,
 *     `ArrangementClosedRightTraits_2`, and
 *     `ArrangementConstractedRightTraits_2`
 */

class ArrangementIdentifiedRightTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_identified_side_tag`.
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

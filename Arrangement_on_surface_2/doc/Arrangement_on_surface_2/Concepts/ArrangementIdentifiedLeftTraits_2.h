/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementIdentifiedLeftTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * identified on the left side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the left boundary side when it is identified.
 *
 * \cgalRefines `ArrangementLeftSideTraits_2`
 *
 * \sa `ArrangementIdentifiedRightTraits_2`,
 *     `ArrangementIdentifiedBottomTraits_2`,
 *     `ArrangementIdentifiedTopTraits_2`,
 *     `ArrangementOpenLeftTraits_2`,
 *     `ArrangementClosedLeftTraits_2`, and
 *     `ArrangementConstractedLeftTraits_2`
 */

class ArrangementIdentifiedLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_identified_side_tag`.
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

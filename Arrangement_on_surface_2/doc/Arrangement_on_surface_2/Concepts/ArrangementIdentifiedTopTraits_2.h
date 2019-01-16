/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementIdentifiedTopTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * identified on the top side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the top boundary side when it is identified.
 *
 * \cgalRefines `ArrangementTopSideTraits_2`
 *
 * \sa `ArrangementIdentifiedLeftTraits_2`,
 *     `ArrangementIdentifiedRightTraits_2`,
 *     `ArrangementIdentifiedBottomTraits_2`,
 *     `ArrangementOpenTopTraits_2`,
 *     `ArrangementClosedTopTraits_2`, and
 *     `ArrangementConstractedTopTraits_2`
 */

class ArrangementIdentifiedTopTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `Arr_identified_side_tag`.
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

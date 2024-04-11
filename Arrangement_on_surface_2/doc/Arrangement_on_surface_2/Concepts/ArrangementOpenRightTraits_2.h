/*!
 * \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementOpenRightTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the right side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * right boundary side when it is open.
 *
 * \cgalRefines{ArrangementRightSideTraits_2}
 *
 * \sa `ArrangementOpenLeftTraits_2`,
 *     `ArrangementOpenBottomTraits_2`,
 *     `ArrangementOpenTopTraits_2`,
 *     `ArrangementClosedRightTraits_2`,
 *     `ArrangementContractedRightTraits_2`, and
 *     `ArrangementIdentifiedVerticalTraits_2`
 */

class ArrangementOpenRightTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
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

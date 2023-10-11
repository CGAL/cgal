/*!
 * \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementOpenTopTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the top side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * top boundary side when it is open.
 *
 * \cgalRefines{ArrangementTopSideTraits_2}
 *
 * \sa `ArrangementOpenLeftTraits_2`,
 *     `ArrangementOpenRightTraits_2`,
 *     `ArrangementOpenBottomTraits_2`,
 *     `ArrangementClosedTopTraits_2`,
 *     `ArrangementContractedTopTraits_2`, and
 *     `ArrangementIdentifiedHorizontalTraits_2`
 */

class ArrangementOpenTopTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
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

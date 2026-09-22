/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosOpenLeftTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the left side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * left boundary side when it is open.
 *
 * \cgalRefines{AosLeftSideTraits_2}
 *
 * \sa `AosOpenRightTraits_2`,
 *     `AosOpenBottomTraits_2`,
 *     `AosOpenTopTraits_2`,
 *     `AosClosedLeftTraits_2`,
 *     `AosContractedLeftTraits_2`, and
 *     `AosIdentifiedVerticalTraits_2`
 */

class AosOpenLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
  typedef unspecified_type Left_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{
  /// @}

  /// \name Accessing Functor Objects
  /// @{
  /// @}
}

/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosOpenRightTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the right side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * right boundary side when it is open.
 *
 * \cgalRefines{AosRightSideTraits_2}
 *
 * \sa `AosOpenLeftTraits_2`,
 *     `AosOpenBottomTraits_2`,
 *     `AosOpenTopTraits_2`,
 *     `AosClosedRightTraits_2`,
 *     `AosContractedRightTraits_2`, and
 *     `AosIdentifiedVerticalTraits_2`
 */

class AosOpenRightTraits_2 {
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

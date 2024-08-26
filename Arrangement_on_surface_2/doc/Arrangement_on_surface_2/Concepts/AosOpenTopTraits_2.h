/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosOpenTopTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the top side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * top boundary side when it is open.
 *
 * \cgalRefines{AosTopSideTraits_2}
 *
 * \sa `AosOpenLeftTraits_2`,
 *     `AosOpenRightTraits_2`,
 *     `AosOpenBottomTraits_2`,
 *     `AosClosedTopTraits_2`,
 *     `AosContractedTopTraits_2`, and
 *     `AosIdentifiedHorizontalTraits_2`
 */

class AosOpenTopTraits_2 {
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

/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosOpenBottomTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is open on
 * the bottom side and curves inserted into the arrangement are expected to
 * reach this boundary side. A model of this concept can handle curves that
 * reach the bottom boundary side when it is open.
 *
 * \cgalRefines{AosBottomSideTraits_2}
 *
 * \sa `AosOpenLeftTraits_2`,
 *     `AosOpenRightTraits_2`,
 *     `AosOpenTopTraits_2`,
 *     `AosClosedBottomTraits_2`,
 *     `AosContractedBottomTraits_2`, and
 *     `AosIdentifiedHorizontalTraits_2`
 */

class AosOpenBottomTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
  typedef unspecified_type Bottom_side_category;
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

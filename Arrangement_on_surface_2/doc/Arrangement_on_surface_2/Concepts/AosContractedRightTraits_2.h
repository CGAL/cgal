/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosContractedRightTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * contracted on the right side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the right boundary side when it is contracted.
 *
 * \cgalRefines{AosRightSideTraits_2}
 *
 * \sa `AosContractedLeftTraits_2`
 * \sa `AosContractedBottomTraits_2`
 * \sa `AosContractedTopTraits_2`
 * \sa `AosOpenRightTraits_2`
 * \sa `AosClosedRightTraits_2`
 * \sa `AosIdentifiedVerticalTraits_2`
 */
class AosContractedRightTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_contracted_side_tag`.
  typedef unspecified_type Right_side_category;
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

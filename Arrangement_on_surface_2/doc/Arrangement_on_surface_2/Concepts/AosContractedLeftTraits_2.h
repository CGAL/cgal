/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosContractedLeftTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * contracted on the left side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the left boundary side when it is contracted.
 *
 * \cgalRefines{AosLeftSideTraits_2}
 *
 * \sa `AosContractedRightTraits_2`
 * \sa `AosContractedBottomTraits_2`
 * \sa `AosContractedTopTraits_2`
 * \sa `AosOpenLeftTraits_2`
 * \sa `AosClosedLeftTraits_2`
 * \sa `AosIdentifiedVerticalTraits_2`
 */
class AosContractedLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_contracted_side_tag`.
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

/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosContractedTopTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is
 * contracted on the top side and curves inserted into the arrangement are
 * expected to reach this boundary side. A model of this concept can handle
 * curves that reach the top boundary side when it is contracted.
 *
 * \cgalRefines{AosTopSideTraits_2}
 *
 * \sa `AosContractedLeftTraits_2`
 * \sa `AosContractedRightTraits_2`
 * \sa `AosContractedBottomTraits_2`
 * \sa `AosOpenTopTraits_2`
 * \sa `AosSlosedTopTraits_2`
 * \sa `AosIdentifiedHorizontalTraits_2`
 */
class AosContractedTopTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_contracted_side_tag`.
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

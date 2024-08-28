/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosClosedBottomTraits_2` must be used when
 * the parameter space of the surface, the arrangement is embedded on, is closed
 * on the left side and curves inserted into the arrangement are expected to
 * reach this boundary side. A model of this concept can handle curves that
 * reach the left boundary side when it is closed.
 *
 * \cgalRefines{AosBottomSideTraits_2}
 *
 * \sa `AosClosedLeftTraits_2`
 * \sa `AosClosedRightTraits_2`
 * \sa `AosClosedTopTraits_2`
 * \sa `AosOpenBottomTraits_2`
 * \sa `AosContractedBottomTraits_2`
 * \sa `AosIdentifiedHorizontalTraits_2`
 */
class AosClosedBottomTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_closed_side_tag`.
  typedef unspecified_type Bottom_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `AosTraits::CompareXOnBoundary_2`.
  typedef unspecified_type Compare_x_on_boundary_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;
  /// @}
}

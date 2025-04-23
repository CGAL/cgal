/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosClosedRightTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is closed on
 * the right side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * right boundary side when it is closed.
 *
 * \cgalRefines{AosRightSideTraits_2}
 *
 * \sa `AosClosedLeftTraits_2`
 * \sa `AosClosedBottomTraits_2`
 * \sa `AosClosedTopTraits_2`
 * \sa `AosOpenRightTraits_2`
 * \sa `AosContractedRightTraits_2`
 * \sa `AosIdentifiedVerticalTraits_2`
 */
class AosClosedRightTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_closed_side_tag`.
  typedef unspecified_type Right_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `AosTraits::CompareYOnBoundary_2`.
  typedef unspecified_type Compare_y_on_boundary_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  /// @}
};

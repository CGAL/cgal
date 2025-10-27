/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosClosedTopTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is closed on
 * the top side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * top boundary side when it is closed.
 *
 * \cgalRefines{AosTopSideTraits_2}
 *
 * \sa `AosClosedLeftTraits_2`
 * \sa `AosClosedRightTraits_2`
 * \sa `AosClosedBottomTraits_2`
 * \sa `AosOpenTopTraits_2`
 * \sa `AosContractedTopTraits_2`
 * \sa `AosIdentifiedHorizontalTraits_2`
 */
class AosClosedTopTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_closed_side_tag`.
  typedef unspecified_type Top_side_category;
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

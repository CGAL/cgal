/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `AosIdentifiedVerticalTraits_2` must be used
 * when the parameter space of the surface, the arrangement is embedded on,
 * is identified on the left and right sides and curves inserted into the
 * arrangement are expected to reach these boundary sides.
 *
 * \cgalRefines{AosBasicTraits_2}
 *
 * \sa `AosIdentifiedHorizontalTraits_2`
 * \sa `AosOpenLeftTraits_2`
 * \sa `AosClosedLeftTraits_2`
 * \sa `AosContractedLeftTraits_2`
 * \sa `AosOpenRightTraits_2`
 * \sa `AosClosedRightTraits_2`
 * \sa `AosContractedRightTraits_2`
 */
class AosIdentifiedVerticalTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_identified_side_tag`.
  typedef unspecified_type Left_side_category;
  typedef unspecified_type Right_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `AosTraits::CompareYOnBoundary_2`.
  typedef unspecified_type Compare_y_on_boundary_2;

  /// models the concept `AosTraits::IsOnYIdentification_2`.
  typedef unspecified_type Is_on_y_identification_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  Is_on_y_identification_2 is_on_y_identification_2_object() const;
  /// @}
}

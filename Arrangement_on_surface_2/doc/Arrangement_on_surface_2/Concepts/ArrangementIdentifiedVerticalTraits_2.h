/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementIdentifiedVerticalTraits_2` must be used
 * when the parameter space of the surface, the arrangement is embedded on, is
 * identified on the left and right sides and curves inserted into the
 * arrangement are expected to reach these boundary sides.
 *
 * \cgalRefines{ArrangementBasicTraits_2}
 *
 * \sa `ArrangementIdentifiedHorizontalTraits_2`,
 *     `ArrangementOpenLeftTraits_2`,
 *     `ArrangementClosedLeftTraits_2`, and
 *     `ArrangementContractedLeftTraits_2`
 *     `ArrangementOpenRightTraits_2`,
 *     `ArrangementClosedRightTraits_2`, and
 *     `ArrangementContractedRightTraits_2`
 */

class ArrangementIdentifiedVerticalTraits_2 {
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

  /// models the concept `ArrTraits::CompareYOnBoundary_2`.
  typedef unspecified_type Compare_y_on_boundary_2;

  /// models the concept `ArrTraits::IsOnYIdentification_2`.
  typedef unspecified_type Is_on_y_identification_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  Is_on_y_identification_2 is_on_y_identification_2_object() const;
  /// @}
}

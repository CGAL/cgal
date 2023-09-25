/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementIdentifiedHorizontalTraits_2` must be used
 * when the parameter space of the surface, the arrangement is embedded on, is
 * identified on the bottom and top sides and curves inserted into the
 * arrangement are expected to reach these boundary sides.
 *
 * \cgalRefines{ArrangementBasicTraits_2}
 *
 * \sa `ArrangementIdentifiedVerticalTraits_2`,
 *     `ArrangementOpenBottomTraits_2`,
 *     `ArrangementClosedBottomTraits_2`, and
 *     `ArrangementContractedBottomTraits_2`
 *     `ArrangementOpenTopTraits_2`,
 *     `ArrangementClosedTopTraits_2`, and
 *     `ArrangementContractedTopTraits_2`
 */

class ArrangementIdentifiedHorizontalTraits_2 {
public:
  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_identified_side_tag`.
  typedef unspecified_type Bottom_side_category;
  typedef unspecified_type Top_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `ArrTraits::CompareXOnBoundary_2`.
  typedef unspecified_type Compare_x_on_boundary_2;

  /// models the concept `ArrTraits::IsOnXIdentification_2`.
  typedef unspecified_type Is_on_x_identification_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Is_on_x_identification_2 is_on_x_identification_2_object() const;
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;
  /// @}
}

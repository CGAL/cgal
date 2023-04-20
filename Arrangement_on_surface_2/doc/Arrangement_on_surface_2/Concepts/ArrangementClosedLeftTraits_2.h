/*!
 * \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * A model of the concept `ArrangementClosedLeftTraits_2` must be used when the
 * parameter space of the surface, the arrangement is embedded on, is closed on
 * the left side and curves inserted into the arrangement are expected to reach
 * this boundary side. A model of this concept can handle curves that reach the
 * left boundary side when it is closed.

 * \cgalRefines{ArrangementLeftSideTraits_2}

\sa `ArrangementClosedRightTraits_2`,
    `ArrangementClosedBottomTraits_2`,
    `ArrangementClosedTopTraits_2`,
    `ArrangementOpenLeftTraits_2`,
    `ArrangementContractedLeftTraits_2`, and
    `ArrangementIdentifiedVerticalTraits_2`,
*/

class ArrangementClosedLeftTraits_2 {
public:
  /// \name Categories
  /// @{

  /*! Must be convertible to `CGAL::Arr_closed_side_tag`.
  */
  typedef unspecified_type Left_side_category;
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `ArrTraits::CompareYOnBoundary_2`.
  typedef unspecified_type Compare_y_on_boundary_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  /// @}
};

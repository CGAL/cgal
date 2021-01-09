namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines AdaptableBinaryFunction
 *
 * \cgalHasModel ArrangementClosedLeftTraits_2::Compare_y_on_boundary_2
 * \cgalHasModel ArrangementClosedRightTraits_2::Compare_y_on_boundary_2
 * \cgalHasModel ArrangementIdentifiedVerticalTraits_2::Compare_y_on_boundary_2
 */
class CompareYOnBoundary_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! Given two points `p1` and `p2` returns `CGAL::SMALLER`, `CGAL::EQUAL`, or
   * `CGAL::LARGER` according to the lexicographic \f$xy\f$-order of the points
   * `p1` and `p2`.
   *
   * \pre \link ArrangementVerticalSideTraits_2::Parameter_space_in_x_2
   * `Parameter_space_in_x_2`\endlink (`p1`) \f$\neq\f$ `CGAL::ARR_INTERIOR` or
   * \link ArrangementVerticalSideTraits_2::Parameter_space_in_x_2
   * `Parameter_space_in_x_2`\endlink (`p2`) \f$\neq\f$ `CGAL::ARR_INTERIOR`.
   */
  Comparison_result operator()(const ArrTraits::Point_2& p1,
                               const ArrTraits::Point_2& p2);

/// @}

}; /* end ArrTraits::CompareYOnBoundary_2 */

}

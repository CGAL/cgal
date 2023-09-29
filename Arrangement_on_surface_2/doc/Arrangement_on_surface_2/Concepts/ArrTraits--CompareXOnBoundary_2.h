namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableFunctor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementClosedBottomTraits_2::Compare_x_on_boundary_2}
 * \cgalHasModels{ArrangementClosedTopTraits_2::Compare_x_on_boundary_2}
 * \cgalHasModels{ArrangementIdentifiedHorizontalTraits_2::Compare_x_on_boundary_2}
 * \cgalHasModelsEnd
 */
class CompareXOnBoundary_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! Given two points `p1` and `p2`, such that either `p1` or `p2` (or both)
   * lie on the bottom or top boundary of the parameter space, compares the \f$
   * x\f$-coordinate of `p1` and the \f$x\f$-coordinate of `p2`. Returns
   * `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` accordingly.
   *
   * \pre \link ArrangementHorizontalSideTraits_2::Parameter_space_in_y_2
   * `Parameter_space_in_y_2`\endlink (`p1`) \f$\neq\f$ `CGAL::ARR_INTERIOR` or
   * \link ArrangementHorizontalSideTraits_2::Parameter_space_in_y_2
   * `Parameter_space_in_y_2`\endlink (`p2`) \f$\neq\f$ `CGAL::ARR_INTERIOR`.
   */
  Comparison_result operator()(const ArrTraits::Point_2& p1,
                               const ArrTraits::Point_2& p2);

/// @}

}; /* end ArrTraits::CompareXOnBoundaryOfCurveEnd_2 */

}

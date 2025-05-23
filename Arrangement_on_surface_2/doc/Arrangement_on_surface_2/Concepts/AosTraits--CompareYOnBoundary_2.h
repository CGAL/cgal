namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableBinaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosClosedLeftTraits_2::Compare_y_on_boundary_2}
 * \cgalHasModels{AosClosedRightTraits_2::Compare_y_on_boundary_2}
 * \cgalHasModels{AosIdentifiedVerticalTraits_2::Compare_y_on_boundary_2}
 * \cgalHasModels{AosSphericalBoundaryTraits_2::Compare_y_on_boundary_2}
 * \cgalHasModelsEnd
 */
class CompareYOnBoundary_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! Given two points `p1` and `p2`, such that either `p1` or `p2` (or both)
   * lie on the bottom or top boundary of the parameter space, compares the
   * \f$y\f$-coordinate of `p1` and the \f$y\f$-coordinate of `p2`. Returns
   * `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` accordingly.
   *
   * \pre \link AosVerticalSideTraits_2::Parameter_space_in_x_2
   * `Parameter_space_in_x_2`\endlink(`p1`) \f$\neq\f$
   * `CGAL::ARR_INTERIOR` or
   * \link AosVerticalSideTraits_2::Parameter_space_in_x_2
   * `Parameter_space_in_x_2`\endlink(`p2`) \f$\neq\f$
   * `CGAL::ARR_INTERIOR`.
   */
  Comparison_result operator()(const AosTraits::Point_2& p1,
                               const AosTraits::Point_2& p2);

  /// @}
}; /* end AosTraits::CompareYOnBoundary_2 */

}

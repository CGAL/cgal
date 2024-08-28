namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableTernaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosOpenBoundaryTraits_2::Compare_x_near_boundary_}
 * \cgalHasModelsEnd
 */
class CompareXNearBoundary_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! Given two \f$x\f$-monotone curves `xcv1` and `xcv2` and an
   * enumeration `ce` that specifies either the minimum ends or the maximum ends
   * of the curves where the curves have a vertical asymptote, compares the
   * \f$x\f$-coordinate of the curves near their respective ends. Returns
   * `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` accordingly. More
   * precisely, compares the \f$x\f$-coordinates of the horizontal projection
   * of a point \f$p\f$ onto `xcv1` and `xcv2`.  If `xcv1` and `xcv2` approach
   * the bottom boundary-side, \f$p\f$ is located far to the bottom, such that
   * the result is invariant under a translation of \f$p\f$ farther to the
   * bottom. If `xcv1` and `xcv2` approach the top boundary-side, \f$p\f$ is
   * located far to the top in a similar manner.
   *
   * \pre The \f$x\f$-coordinates of the boundary of the curves at their
   * respective ends are equal. That is,
   * \link AosHorizontalSideTraits_2::Compare_x_on_boundary_2
   * `compare_x_on_boundary_2`\endlink(`xcv1`, `xcv2`, `ce`) = `CGAL::EQUAL`.
   *
   * \pre \link AosHorizontalSideTraits_2::Parameter_space_in_y_2
   * `parameter_space_in_y_2`\endlink(`xcv1`, `ce`) =
   * \link AosHorizontalSideTraits_2::Parameter_space_in_y_2
   * `parameter_space_in_y_2`\endlink(`xcv2`, `ce`).
   *
   * \pre \link AosHorizontalSideTraits_2::Parameter_space_in_y_2
   * `parameter_space_in_y_2`\endlink(`xcv1`, `ce`) \f$\neq\f$
   * `CGAL::ARR_INTERIOR`.
   */
  Comparison_result operator()(const AosTraits::X_monotone_curve_2& xcv1,
                               const AosTraits::X_monotone_curve_2& xcv2,
                               CGAL::Arr_curve_end ce);

  /// @}
}; /* end AosTraits::CompareXNearBoundary_2 */

}

namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableFunctor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementHorizontalSideTraits_2::Compare_x_on_boundary_2}
 * \cgalHasModels{ArrangementOpenBoundaryTraits_2::Compare_x_on_boundary_2}
 * \cgalHasModels{ArrangementSphericalBoundaryTraits_2::Compare_x_on_boundary_2}
 * \cgalHasModelsEnd
 */
class CompareXOnBoundaryOfCurveEnd_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! Given a point `p`, an \f$x\f$-monotone curve `xcv`, and an
   * enumeration `ce` that specifies either the minimum or the maximum end of
   * the curve where the curve has a vertical asymptote, compares the \f$
   * x\f$-coordinate of `p` and the \f$x\f$-coordinate of the limit of the
   * curve at its specified end. The variable `xcv` identifies the parametric
   * curve \f$c(t) = (x(t), y(t))\f$ defined over an open or half-open interval
   * with endpoints \f$ 0\f$ and \f$ 1\f$. The enumeration `ce` identifies an
   * open end \f$d \in\{0,1\}\f$ of \f$c\f$.  Formally, compares the \f$
   * x\f$-coordinate of `p` and \f$ \lim_{t \rightarrow d} x(t)\f$. Returns
   * `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` accordingly.
   *
   * \pre \link ArrangementHorizontalSideTraits_2::Parameter_space_in_y_2 `Parameter_space_in_y_2`\endlink (`xcv`, `ce`) \f$\neq\f$ `CGAL::ARR_INTERIOR`.
   *
   * \pre If the parameter space is unbounded, \f$c\f$ has a vertical asymptote
   *      at its \f$ d\f$-end; that is,
   *      \link ArrangementVerticalSideTraits_2::Parameter_space_in_x_2 `Parameter_space_in_x_2`\endlink(`xcv`, `ce`) = `CGAL::ARR_INTERIOR`.
   */
  Comparison_result operator()(const ArrTraits::Point_2& p,
                               const ArrTraits::X_monotone_curve_2& xcv,
                               CGAL::Arr_curve_end ce);

/*! Given two \f$ x\f$-monotone curves `xcv1` and `xcv2` and two indices `ce1`
 * and `ce2` that specify either the minimum or the maximum ends of `xcv1` and
 * `xcv2`, respectively, where the curves have vertical asymptotes, compares the
 * \f$ x\f$-coordinates of the limits of the curves at their specified
 * ends. The variables `xcv1` and `xcv2` identify the parametric curves \f$
 * c_1(t) = (x_1(t),y_1(t))\f$ and \f$ c_2(t) = (x_2(t),y_2(t))\f$,
 * respectively, defined over open or half-open intervals with endpoints \f$
 * 0\f$ and \f$ 1\f$. The indices `ce1` and `ce2` identify open ends \f$ d_1
 * \in\{0,1\}\f$ and \f$ d_2 \in\{0,1\}\f$ of \f$ c_1\f$ and \f$ c_2\f$,
 * respectively. Formally, compares \f$ \lim_{t \rightarrow d_1} x_1(t)\f$ and
 * \f$\lim_{t \rightarrow d_2} x_2(t)\f$. Returns `CGAL::SMALLER`,
 * `CGAL::EQUAL`, or `CGAL::LARGER` accordingly.
 *
 * \pre \link ArrangementHorizontalSideTraits_2::Parameter_space_in_y_2 `Parameter_space_in_y_2`\endlink(`xcv1`, `ce1`) \f$\neq\f$ `CGAL::ARR_INTERIOR`.
 *
 * \pre \link ArrangementHorizontalSideTraits_2::Parameter_space_in_y_2 `Parameter_space_in_y_2`\endlink(`xcv2`, `ce2`) \f$\neq\f$ `CGAL::ARR_INTERIOR`.
 *
 * \pre If the parameter space is unbounded, \f$c_1\f$ has a vertical
 *      asymptote at its respective end; that is,
 *      \link ArrangementVerticalSideTraits_2::Parameter_space_in_x_2 `Parameter_space_in_x_2`\endlink(`xcv1`, `ce1`) = `CGAL::ARR_INTERIOR`.
 *
 * \pre If the parameter space is unbounded, \f$c_2\f$ has a vertical asymptote
 *      at its respective end; that is,
 *      \link ArrangementVerticalSideTraits_2::Parameter_space_in_x_2 `Parameter_space_in_x_2`\endlink(`xcv2`, `ce2`) = `CGAL::ARR_INTERIOR`.
 */
Comparison_result operator()(const ArrTraits::X_monotone_curve_2& xcv1,
                             CGAL::Arr_curve_end ce1,
                             const ArrTraits::X_monotone_curve_2& xcv2,
                             CGAL::Arr_curve_end ce2);

/// @}

}; /* end ArrTraits::CompareXOnBoundaryOfCurveEnd_2 */

}

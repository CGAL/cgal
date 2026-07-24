namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableTernaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosBasicTraits_2::Compare_y_at_x_right_2}
 * \cgalHasModelsEnd
 */
class CompareYAtXRight_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! accepts two \f$x\f$-monotone curves `xc1` and `xc2`
   * that have a common left endpoint `p`, and returns `CGAL::SMALLER`,
   * `CGAL::EQUAL` or `CGAL::LARGER` according to the relative position of the
   * two curves immediately to the right of \f$p\f$. Note that in case one of
   * the \f$x\f$-monotone curves is a vertical segment emanating upward from
   * `p`, it is always considered to be <I>above</I> the other curve.
   */
  Comparison_result operator()(AosTraits::X_monotone_curve_2 xc1,
                               AosTraits::X_monotone_curve_2 xc2,
                               AosTraits::Point_2 p);

  /// @}
}; /* end AosTraits::CompareYAtXRight_2 */

}

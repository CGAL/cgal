namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableBinaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosBasicTraits_2::Compare_y_at_x_2}
 * \cgalHasModelsEnd
 */
class CompareYAtX_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! compares the \f$y\f$-coordinates of `p` and the vertical projection
   * of `p` on `xc`, and returns `CGAL::SMALLER`, `CGAL::EQUAL`, or
   * `CGAL::LARGER` according to the result.
   */
  Comparison_result operator()(AosTraits::Point_2 p,
                               AosTraits::X_monotone_curve_2 xc);

  /// @}
}; /* end AosTraits::CompareYAtX_2 */

}

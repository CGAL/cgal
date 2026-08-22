namespace AocTraits {

/*! \ingroup PkgArrangementOnCurve1ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AocTraits_1::Compare_x_1}
 * \cgalHasModelsEnd
 */
class CompareX_1 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! returns `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` according to the
   * \f$x\f$-ordering of points `p1` and `p2`.
   */
  Comparison_result operator()(AosTraits::Point_1 p1, AosTraits::Point_1 p2);

  /// @}
};

}

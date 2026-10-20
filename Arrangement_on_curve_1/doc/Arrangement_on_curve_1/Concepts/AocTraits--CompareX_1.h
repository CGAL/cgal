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

  /*! returns `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` according to
   * the ordering of points `p1` and `p2` along the 1D carrier curve.
   */
  Comparison_result operator()(const AocTraits::Point_1& p1, const AocTraits::Point_1& p2);

  /// @}
};

}

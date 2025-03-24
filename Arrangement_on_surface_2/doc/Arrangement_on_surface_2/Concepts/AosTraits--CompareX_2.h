namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableBinaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosBasicTraits_2::Compare_x_2}
 * \cgalHasModelsEnd
 */
class CompareX_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! returns `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` according to the
   * \f$x\f$-ordering of points `p1` and `p2`.
   */
  Comparison_result operator()(AosTraits::Point_2 p1, AosTraits::Point_2 p2);

  /// @}
}; /* end AosTraits::CompareX_2 */

}

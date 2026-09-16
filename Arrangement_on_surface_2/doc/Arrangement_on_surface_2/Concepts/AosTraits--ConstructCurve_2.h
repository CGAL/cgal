namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosConstructCurveTraits_2::Construct_curve_2}
 * \cgalHasModelsEnd
 */
class ConstructCurve_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! returns a curve connecting `p1` and `p2` (i.e., the two input points
   * are its endpoints).
   */
  AosTraits::Curve_2 operator()(AosTraits::Point_2 p1, AosTraits::Point_2 p2);
  /// @}
}; /* end AosTraits::ConstructCurve_2 */

}

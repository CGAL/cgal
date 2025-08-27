namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosIdentifiedVerticalTraits_2::Is_on_y_identification_2}
 * \cgalHasModelsEnd
 */
class IsOnYIdentification_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines whether the point `p` lies on the horizontal identification
   * curve.
   */
  bool operator()(AosTraits::Point& point) const;

  /*! determines whether the curve `cv` lies on the horizontal identification
   * curve.
   */
  bool operator()(AosTraits::const X_monotone_curve_2& cv) const;

  /// @}
}; /* end AosTraits::IsOnYIdentification_2 */

}

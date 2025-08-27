namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosIdentifiedHorizontalTraits_2::Is_on_x_identification_2}
 * \cgalHasModelsEnd
 */
class IsOnXIdentification_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines whether the point `p` lies on the vertical identification
   * curve.
   */
  bool operator()(AosTraits::Point& point) const;

  /*! determines whether the curve `cv` lies on the vertical identification
   * curve.
   */
  bool operator()(AosTraits::const X_monotone_curve_2& cv) const;

  /// @}
}; /* end AosTraits::IsOnXIdentification_2 */

}

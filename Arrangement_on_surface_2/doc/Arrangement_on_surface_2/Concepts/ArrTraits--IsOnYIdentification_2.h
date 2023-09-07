namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementIdentifiedVerticalTraits_2::Is_on_y_identification_2}
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
  bool operator()(ArrTraits::Point& point) const;

  /*! determines whether the curve `cv` lies on the horizontal identification
   * curve.
   */
  bool operator()(ArrTraits::const X_monotone_curve_2& cv) const;

  /// @}

}; /* end ArrTraits::IsOnYIdentification_2 */

}

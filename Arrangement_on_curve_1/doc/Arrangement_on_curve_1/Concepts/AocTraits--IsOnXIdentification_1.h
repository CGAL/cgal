namespace AocTraits {

/*! \ingroup PkgArrangementOnCurve1ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosIdentifiedTraits_1::Is_on_x_identification_1}
 * \cgalHasModelsEnd
 */
class IsOnXIdentification_1 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines whether the point `p` coincides with the identification point.
   */
  bool operator()(AosTraits::Point& point) const;

  /*! determines whether the point `p` coincides with the identification curvepoint.
   */
  bool operator()(AosTraits::const Point_1& p) const;

  /// @}
}; /* end AosTraits::IsOnXIdentification_1 */

}

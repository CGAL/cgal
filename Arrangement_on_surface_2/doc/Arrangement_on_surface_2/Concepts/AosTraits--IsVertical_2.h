namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosBasicTraits_2::Is_vertical_2}
 * \cgalHasModelsEnd
 */
class IsVertical_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines whether `xc` is a vertical segment.
   */
  bool operator()(AosTraits::X_monotone_curve_2 xc);

  /// @}
}; /* end AosTraits::IsVertical_2 */

}

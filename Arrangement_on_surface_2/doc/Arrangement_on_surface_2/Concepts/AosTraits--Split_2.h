namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosXMonotoneTraits_2::Split_2}
 * \cgalHasModelsEnd
 */
class Split_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! accepts an input curve `xc` and a split point `p` in its interior. It
   * splits `xc` at the split point into two subcurves `xc1` and `xc2`, such
   * that `p` is `xc1`'s <I>right</I> endpoint and `xc2`'s <I>left</I> endpoint.
   */
  void operator()(AosTraits::X_monotone_curve_2 xc,
                  AosTraits::Point_2 p,
                  AosTraits::X_monotone_curve_2& xc1,
                  AosTraits::X_monotone_curve_2& xc2);

  /// @}
}; /* end AosTraits::Split_2 */

}

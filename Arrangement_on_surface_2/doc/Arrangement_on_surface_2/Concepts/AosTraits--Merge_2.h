namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosXMonotoneTraits_2::Merge_2}
 * \cgalHasModelsEnd
 */
class Merge_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! accepts two <I>mergeable</I> \f$x\f$-monotone curves `xc1` and `xc2`
   * and assigns `xc` with the merged curve.
   *
   * \pre `are_mergeable_2`(`xc1`, `xc2`) is true.
   */
  void merge(AosTraits::X_monotone_curve_2 xc1,
             AosTraits::X_monotone_curve_2 xc2,
             AosTraits::X_monotone_curve_2& xc);

  /// @}
}; /* end AosTraits::Merge_2 */

}

namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosApproximatePointTraits_2::Approximate_2}
 * \cgalHasModels{AosApproximateTraits_2::Approximate_2}
 * \cgalHasModelsEnd
 */
class ApproximatePoint_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! obtains an approximation of `p`'s \f$x\f$-coordinate (if `i == 0`), or of
   * `p`'s \f$y\f$-coordinate (if `i == 1`).
   * \pre `i` is either 0 or 1.
   */
  Approximate_number_type operator()(AosTraits::Point_2 p, int i);

  /// @}
}; /* end AosTraits::Approximate_2 */

}

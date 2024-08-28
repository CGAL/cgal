namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosApproximatePointTraits_2::Approximate_2}
 * \cgalHasModelsEnd
 */
class Approximate_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! obtains an approximation of `p`'s \f$x\f$-coordinate (if `i == 0`), or of
   * `p`'s \f$y\f$-coordinate (if `i == 1`).
   */
  CGAL::Approximate_number_type operator()(AosTraits::Point_2 p, int i);

  /// @}
}; /* end AosTraits::Approximate_2 */

}

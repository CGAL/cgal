namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementXMonotoneTraits_2::Are_mergeable_2}
 * \cgalHasModelsEnd
 */
class AreMergeable_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! accepts two \f$x\f$-monotone curves `xc1` and `xc2` and determines
   * whether they can be merged to form a single \f$x\f$-monotone curve. `xc1`
   * and `xc2` are mergeable if their underlying curves are identical, they
   * share a common endpoint, and they do not bend to form a
   * non-\f$x\f$-monotone curve.
   */
  bool operator()(ArrTraits::X_monotone_curve_2 xc1,
                  ArrTraits::X_monotone_curve_2 xc2);

  /// @}

}; /* end ArrTraits::AreMergeable_2 */

}

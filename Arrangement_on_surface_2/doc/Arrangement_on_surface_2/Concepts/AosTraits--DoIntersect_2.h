namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosXMonotoneTraits_2::DoIntersect_2}
 * \cgalHasModelsEnd
 */
class DoIntersect_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines whether two \f$x\f$-monotone curves intersect.
   *
   * \param xc1 The first \f$x\f$-monotone curve.
   * \param xc2 The second \f$x\f$-monotone curve.
   * \param closed a boolean flag that indicates whether `xc1` and `xc2` are closed; the defalut is `true`.
   * \return `true` if `xc1` and `xc2` are closed and intersect or if `xc1` and `xc2` are open and intersect
   * in at least one of their interiors, and `false` otherwise.
   */
  template <typename OutputIterator>
  OutputIterator operator()(AosTraits::X_monotone_curve_2 xc1,
                            AosTraits::X_monotone_curve_2 xc2,
                            bool closed = true);

  /// @}

}; /* end AosTraits::DoIntersect_2 */

}

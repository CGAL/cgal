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
   * \param consider_common_endpoints indicates whether common endpoints should be counted as intersections.
   * \return `true` if `consider_common_endpoints` is `true` and `xcv1` and `xcv2` intersect or if
   *  `consider_common_endpoints` is `false` and at least one of the interiors of `xcv1` and `xcv2` intersect,
   *   and `false` otherwise.
   * in at least one of their interiors, and `false` otherwise.
   */
  template <typename OutputIterator>
  OutputIterator operator()(AosTraits::X_monotone_curve_2 xc1,
                            AosTraits::X_monotone_curve_2 xc2,
                            bool consider_common_endpoints = true);

  /// @}

}; /* end AosTraits::DoIntersect_2 */

}

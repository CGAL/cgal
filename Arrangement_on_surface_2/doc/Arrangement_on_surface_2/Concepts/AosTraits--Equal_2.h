namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableBinaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosBasicTraits_2::Equal_2}
 * \cgalHasModelsEnd
 */
class Equal_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines whether `p1` and `p2` are geometrically equivalent.
   */
  bool operator()(AosTraits::Point_2 p1, AosTraits::Point_2 p2);

  /*! determines whether `xc1` and `xc2` are geometrically equivalent
   * (have the same graph).
   */
  bool operator()(AosTraits::X_monotone_curve_2 xc1,
                  AosTraits::X_monotone_curve_2 xc2);

  /// @}
}; /* end AosTraits::Equal_2 */

}

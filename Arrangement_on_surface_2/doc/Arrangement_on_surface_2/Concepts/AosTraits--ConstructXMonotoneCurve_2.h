namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosConstructXMonotoneCurveTraits_2::Construct_x_monotone_curve_2}
 * \cgalHasModelsEnd
 */
class ConstructXMonotoneCurve_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! returns an \f$x\f$-monotone curve connecting `p1` and `p2` (i.e., the
   * two input points are its endpoints).
   */
  AosTraits::X_monotone_curve_2 operator()(AosTraits::Point_2 p1,
                                           AosTraits::Point_2 p2);

  /// @}
}; /* end AosTraits::ConstructXMonotoneCurve_2 */

}

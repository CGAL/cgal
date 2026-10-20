namespace AocTraits {

/*! \ingroup PkgArrangementOnCurve1ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AocClosedBoundaryTraits_1::Parameter_space_in_x_1}
 * \cgalHasModelsEnd
 */
class ParameterSpaceInX_1 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! determines the location of a point pre-image \f$t\f$ along the interval
   * parameter space.
   * If \f$t\f$ lies on the left boundary, returns `CGAL::Arrangement_on_curve_1::LEFT_BOUNDARY`;
   * otherwise, if \f$t\f$ lies on the right boundary, returns `CGAL::Arrangement_on_curve_1::RIGHT_BOUNDARY`;
   * otherwise, returns `CGAL::Arrangement_on_curve_1::INTERIOR`.
   */
  Parameter_space operator()(const AocTraits::Point_1& p);

  /// @}
};

}

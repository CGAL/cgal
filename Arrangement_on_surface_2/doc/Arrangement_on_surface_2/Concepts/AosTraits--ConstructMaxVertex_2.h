namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosBasicTraits_2::Construct_max_vertex_2}
 * \cgalHasModelsEnd
 */
class ConstructMaxVertex_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! returns the lexicographically larger (right) endpoint of `xc`.
   */
  AosTraits::Point_2 operator()(AosTraits::X_monotone_curve_2 xc);

  /// @}
}; /* end AosTraits::ConstructMaxVertex_2 */

}

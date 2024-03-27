namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{AdaptableUnaryFunction}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementBasicTraits_2::Construct_min_vertex_2}
 * \cgalHasModelsEnd
 */
class ConstructMinVertex_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! returns the lexicographically smaller (left) endpoint of `xc`.
   */
  ArrTraits::Point_2 operator()(ArrTraits::X_monotone_curve_2 xc);

  /// @}

}; /* end ArrTraits::ConstructMinVertex_2 */

}

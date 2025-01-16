namespace ArrDirectionalTraits {

/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines{AdaptableUnaryFunction}

\cgalHasModelsBegin
\cgalHasModels{ArrangementDirectionalXMonotoneTraits_2::Construct_opposite_2}
\cgalHasModelsEnd

*/

class ConstructOpposite_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts an \f$ x\f$-monotone curve `xc` and returns its opposite curve,
namely a curve whose graph is the same as `xc`'s, and whose source and
target are swapped with respect to `xc`'s source and target.
*/
ArrDirectionalTraits::X_monotone_curve_2 operator()(ArrDirectionalTraits::X_monotone_curve_2 xc);

/// @}

}; /* end ArrDirectionalTraits::ConstructOpposite_2 */

}

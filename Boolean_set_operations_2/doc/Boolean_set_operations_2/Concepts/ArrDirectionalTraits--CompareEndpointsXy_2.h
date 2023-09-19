namespace ArrDirectionalTraits {

/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines{AdaptableUnaryFunction}

\cgalHasModelsBegin
\cgalHasModels{ArrangementDirectionalXMonotoneTraits_2::Compare_endpoints_xy_2}
\cgalHasModelsEnd

*/

class CompareEndpointsXy_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts an input curve `xc` and compares its source and target
points. It returns `SMALLER` if the curve is directed from
lexicographically left to right, and `LARGER` if it is directed
from lexicographically right to left.
*/
Comparison_result operator()(ArrDirectionalTraits::X_monotone_curve_2 xc);

/// @}

}; /* end ArrDirectionalTraits::CompareEndpointsXy_2 */

}

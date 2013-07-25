namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableTernaryFunction 

\cgalHasModel ArrangementOpenBoundaryTraits_2::Compare_y_near_boundary_2 

*/

class CompareYNearBoundary_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Given two \f$ x\f$-monotone curves `xcv1` and `xcv2` and an 
enumeration `ce` that specifies either the minimum or the maximum 
ends of the curves, compares the \f$ y\f$-coordinate of the curves near 
their respective ends. Returns `SMALLER`, `EQUAL`, or 
`LARGER` accordingly. More precisely, compares the 
\f$ y\f$-coordinates of the vertical projection of a point \f$ p\f$ onto 
predicate `Parameter_space_in_x_2` evaluates to 
`ARR_LEFT_BOUNDARY` when applied to `xcv1` and `ce` and 
when applied to `xcv2` and `ce`. In this case \f$ p\f$ is 
located far to the left, such that the result is invariant under 
a translation of \f$ p\f$ farther to the left. If `ce` is 
evaluates to `ARR_RIGHT_BOUNDARY` when applied to `xcv1` and 
`ce` and when applied to `xcv2` and `ce`. In that case 
\f$ p\f$ is located far to the right in a similar manner. 
\pre `parameter_space_in_x_2`(`xcv2`, `ce`) = `parameter_space_in_x_2`(`xcv1`, `ce`). 
\pre `parameter_space_in_x_2`(`xcv1`, `ce`) \f$ \neq\f$ `ARR_INTERIOR`. 
*/ 
Comparison_result operator()(const ArrTraits::X_monotone_curve_2& xcv1, 
const ArrTraits::X_monotone_curve_2& xcv2, 
Arr_curve_end ce); 

/// @}

}; /* end ArrTraits::CompareYNearBoundary_2 */

}

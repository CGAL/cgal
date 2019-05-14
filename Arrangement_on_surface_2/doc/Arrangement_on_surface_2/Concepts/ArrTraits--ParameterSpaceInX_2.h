namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableBinaryFunction 

\cgalHasModel ArrangementOpenBoundaryTraits_2::Parameter_space_in_x_2 

*/

class ParameterSpaceInX_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!

Given an \f$ x\f$-monotone curve `xcv` and an enumeration `ce` 
that specifies either the minimum or the maximum end of the curve, 
determines the location of the curve end along the \f$ x\f$-dimension. 
The variable `xcv` identifies the parametric curve 
\f$ C(t) = (X(t),Y(t))\f$ defined over an open or half-open interval with 
endpoints \f$ 0\f$ and \f$ 1\f$. The enumeration `ce` identifies an open 
end \f$ d \in\{0,1\}\f$ of \f$ C\f$. Formally, determines whether 
\f$ \lim_{t \rightarrow d} X(t)\f$ evaluates to \f$ b_l\f$, \f$ b_r\f$, or a value 
in between, where \f$ b_l\f$ and \f$ b_r\f$ are the \f$ x\f$-coordinates of the 
left and right boundaries of the parameter space, respectively. 
Returns `ARR_LEFT_BOUNDARY`, `ARR_RIGHT_BOUNDARY`, or 
`ARR_INTERIOR`, accordingly. 
*/ 
Arr_parameter_space operator()(const ArrTraits::X_monotone_curve_2& xcv, 
Arr_curve_end ce); 

/// @}

}; /* end ArrTraits::ParameterSpaceInX_2 */

}

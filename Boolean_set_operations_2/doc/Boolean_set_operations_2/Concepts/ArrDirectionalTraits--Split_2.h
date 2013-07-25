namespace ArrDirectionalTraits {
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines `AdaptableUnaryFunction`

\cgalHasModel `ArrangementDirectionalXMonotoneTraits_2::Split_2` 

*/

class Split_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts an input curve `xc` and a split point `p` in its 
interior. It splits `xc` at the split point into two subcurves 
`xc1` and `xc2`, such that `p` is `xc1`'s <I>right</I> 
endpoint and `xc2`'s <I>left</I> endpoint. The direction of `xc` 
is preserved. That is, in case `xc` is directed from left to right, 
`p` becomes `xc1`'s target and `c2`'s source; 
otherwise, `p` becomes `xc2`'s target and `xc1`'s source. 
*/ 
void operator()(ArrDirectionalTraits::X_monotone_curve_2 xc, 
ArrDirectionalTraits::Point_2 p, 
ArrDirectionalTraits::X_monotone_curve_2& xc1, 
ArrDirectionalTraits::X_monotone_curve_2& xc2); 

/// @}

}; /* end ArrDirectionalTraits::Split_2 */

}
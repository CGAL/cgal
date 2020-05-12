namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor

\cgalHasModel ArrangementXMonotoneTraits_2::Split_2

*/

class Split_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts an input curve `xc` and a split point `p` in its
interior. It splits `xc` at the split point into two subcurves `xc1`
and `xc2`, such that `p` is `xc1`'s <I>right</I> endpoint and
`xc2`'s <I>left</I> endpoint.
*/
void operator()(ArrTraits::X_monotone_curve_2 xc,
ArrTraits::Point_2 p,
ArrTraits::X_monotone_curve_2& xc1,
ArrTraits::X_monotone_curve_2& xc2);

/// @}

}; /* end ArrTraits::Split_2 */

}

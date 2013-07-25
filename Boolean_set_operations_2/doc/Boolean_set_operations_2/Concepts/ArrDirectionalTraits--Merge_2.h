namespace ArrDirectionalTraits {
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines `AdaptableBinaryFunction` 

\cgalHasModel `ArrangementDirectionalXMonotoneTraits_2::Merge_2` 

*/

class Merge_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts two <I>mergeable</I> \f$ x\f$-monotone curves `xc1` and 
`xc2` and asigns `xc` with the merged curve. If the target 
point of `xc1` and the source point of `xc2` coincide; then 
the source point of `xc1` and the target point of `xc2` become 
the source and target points of `xc`, respectively. If the target 
point of `xc2` and the source point of `xc1` coincide; then 
the source point of `xc2` and the target point of `xc1` become 
the source and target points of `xc`, respectively. 
\pre `are_mergeable_2`(`xc1`, `xc2`) is true. 
*/ 
void operator()(ArrDirectionalTraits::X_monotone_curve_2 xc1, 
ArrDirectionalTraits::X_monotone_curve_2 xc2, 
ArrDirectionalTraits::X_monotone_curve_2& xc); 

/// @}

}; /* end ArrDirectionalTraits::Merge_2 */

}

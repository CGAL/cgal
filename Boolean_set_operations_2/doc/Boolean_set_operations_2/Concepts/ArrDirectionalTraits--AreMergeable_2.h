namespace ArrDirectionalTraits {
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines `AdaptableBinaryFunction` 

\cgalHasModel `ArrangementDirectionalXMonotoneTraits_2::Are_mergeable_2` 

*/

class AreMergeable_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts two \f$ x\f$-monotone curves `xc1` and `xc2` and determines 
whether they can be merged to form a single \f$ x\f$-monotone curve. 
`xc1` and `xc2` are mergeable if their underlying curves are 
identical, they share a common endpoint, and they do not bend to form 
a non-\f$ x\f$-monotone curve. 
\pre The target point of `xc1` and the source point `xc2` coincide or the source point of `xc2` and the target point `xc2` coincide. 
*/ 
bool operator()(ArrDirectionalTraits::X_monotone_curve_2 xc1, 
ArrDirectionalTraits::X_monotone_curve_2 xc2); 

/// @}

}; /* end ArrDirectionalTraits::AreMergeable_2 */

}
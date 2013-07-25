namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor 

\cgalHasModel ArrangementXMonotoneTraits_2::Merge_2 

*/

class Merge_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
accepts two <I>mergeable</I> \f$ x\f$-monotone curves `xc1` and `xc2` 
and asigns `xc` with the merged curve. 
\pre `are_mergeable_2`(`xc1`, `xc2`) is true. 
*/ 
void merge(ArrTraits::X_monotone_curve_2 xc1, 
ArrTraits::X_monotone_curve_2 xc2, 
ArrTraits::X_monotone_curve_2& xc); 

/// @}

}; /* end ArrTraits::Merge_2 */

}

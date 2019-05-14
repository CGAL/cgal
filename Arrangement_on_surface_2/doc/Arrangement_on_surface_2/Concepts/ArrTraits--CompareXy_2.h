namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableBinaryFunction

\cgalHasModel ArrangementBasicTraits_2::Compare_xy_2 

*/

class CompareXy_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns `SMALLER`, `EQUAL`, or `LARGER` according 
to the lexicographic \f$ xy\f$-order of the points `p1` and `p2`. 
*/ 
Comparison_result operator()(ArrTraits::Point_2 p1, 
ArrTraits::Point_2 p2); 

/// @}

}; /* end ArrTraits::CompareXy_2 */

}

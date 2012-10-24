namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalconcept

\refines AdaptableBinaryFunction 

\hasModel ArrangementBasicTraits_2::Compare_x_2 

*/

class CompareX_2 {
public:

/// \name Has Models 
/// @{

/*! 
returns `SMALLER`, `EQUAL`, or `LARGER` 
according to the \f$ x\f$-ordering of points `p1` and `p2`. 
*/ 
Comparison_result operator()(ArrTraits::Point_2 p1, 
ArrTraits::Point_2 p2); 

/// @}

}; /* end ArrTraits::CompareX_2 */

}

namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalconcept

\refines AdaptableUnaryFunction 

\hasModel ArrangementBasicTraits_2::Is_vertical_2 

*/

class IsVertical_2 {
public:

/// \name Has Models 
/// @{

/*! 
determines whether `xc` is a vertical segment. 
*/ 
bool operator()(ArrTraits::X_monotone_curve_2 xc); 

/// @}

}; /* end ArrTraits::IsVertical_2 */

}

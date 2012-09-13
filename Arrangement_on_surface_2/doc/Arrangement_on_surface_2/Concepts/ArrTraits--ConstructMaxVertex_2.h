
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalconcept

\refines AdaptableUnaryFunction 

\hasModel ArrangementBasicTraits_2::Construct_max_vertex_2 

*/

class ArrTraits::ConstructMaxVertex_2 {
public:

/// \name Has Models 
/// @{

/*! 
returns the lexicographically larger (right) endpoint of `xc`. 
*/ 
ArrTraits::Point_2 operator()(ArrTraits::X_monotone_curve_2 xc); 

/// @}

}; /* end ArrTraits::ConstructMaxVertex_2 */


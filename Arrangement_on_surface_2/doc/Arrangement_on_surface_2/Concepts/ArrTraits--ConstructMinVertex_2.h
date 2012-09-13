
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalconcept

\refines AdaptableUnaryFunction 

\hasModel ArrangementBasicTraits_2::Construct_min_vertex_2 

*/

class ArrTraits::ConstructMinVertex_2 {
public:

/// \name Has Models 
/// @{

/*! 
returns the lexicographically smaller (left) endpoint of `xc`. 
*/ 
ArrTraits::Point_2 operator()(ArrTraits::X_monotone_curve_2 xc); 

/// @}

}; /* end ArrTraits::ConstructMinVertex_2 */


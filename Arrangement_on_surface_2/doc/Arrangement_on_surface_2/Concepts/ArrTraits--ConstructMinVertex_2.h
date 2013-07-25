namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableUnaryFunction 

\cgalHasModel ArrangementBasicTraits_2::Construct_min_vertex_2 

*/

class ConstructMinVertex_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the lexicographically smaller (left) endpoint of `xc`. 
*/ 
ArrTraits::Point_2 operator()(ArrTraits::X_monotone_curve_2 xc); 

/// @}

}; /* end ArrTraits::ConstructMinVertex_2 */

}

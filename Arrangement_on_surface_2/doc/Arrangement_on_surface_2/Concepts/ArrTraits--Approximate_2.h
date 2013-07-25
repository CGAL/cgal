namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor 

\cgalHasModel ArrangementLandmarkTraits_2::Approximate_2 

*/

class Approximate_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns an approximation of `p`'s \f$ x\f$-coordinate (if `i == 0`), 
or of `p`'s \f$ y\f$-coordinate (if `i == 1`). 
*/ 
Approximate_number_type operator()( ArrTraits::Point_2 p, 
int i); 

/// @}

}; /* end ArrTraits::Approximate_2 */

}

namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableUnaryFunction

\cgalHasModel ArrangementBasicTraits_2::Construct_max_vertex_2

*/

class ConstructMaxVertex_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the lexicographically larger (right) endpoint of `xc`.
*/
ArrTraits::Point_2 operator()(ArrTraits::X_monotone_curve_2 xc);

/// @}

}; /* end ArrTraits::ConstructMaxVertex_2 */

}

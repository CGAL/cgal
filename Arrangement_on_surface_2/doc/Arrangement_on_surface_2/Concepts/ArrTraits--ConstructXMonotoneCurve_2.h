namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor

\cgalHasModel ArrangementLandmarkTraits_2::Construct_x_monotone_curve_2

*/

class ConstructXMonotoneCurve_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns an \f$ x\f$-monotone curve connecting `p1` and `p2` (i.e., the
two input points are its endpoints).
*/
ArrTraits::X_monotone_curve_2 operator() ( ArrTraits::Point_2 p1,
ArrTraits::Point_2 p2);

/// @}

}; /* end ArrTraits::ConstructXMonotoneCurve_2 */

}

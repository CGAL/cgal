namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableBinaryFunction

\cgalHasModel ArrangementBasicTraits_2::Compare_x_2

*/

class CompareX_2 {
public:

/// \name Operations
/// A model of this concept must provide:
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

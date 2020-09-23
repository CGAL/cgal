namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines AdaptableBinaryFunction

\cgalHasModel ArrangementBasicTraits_2::Compare_y_at_x_2

*/

class CompareYAtX_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
compares the \f$ y\f$-coordinates of `p` and the vertical
projection of `p` on `xc`, and returns `SMALLER`, `EQUAL`,
or `LARGER` according to the result.
*/
Comparison_result operator()(ArrTraits::Point_2 p,
ArrTraits::X_monotone_curve_2 xc);

/// @}

}; /* end ArrTraits::CompareYAtX_2 */

}


/*!
\ingroup PkgArrangement2Concepts
\cgalconcept

\refines AdaptableBinaryFunction

\hasModel ArrangementBasicTraits_2::Compare_xy_2 

*/

class ArrTraits::CompareXy_2 {
public:

/// \name Has Models 
/// @{

/*! 
returns `SMALLER`, `EQUAL`, or `LARGER` according 
to the lexicographic \f$ xy\f$-order of the points `p1` and `p2`. 
*/ 
Comparison_result operator()(ArrTraits::Point_2 p1, 
ArrTraits::Point_2 p2); 

/// @}

}; /* end ArrTraits::CompareXy_2 */



/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalconcept

\refines ::AdaptableUnaryFunctor 

\hasModel `ArrangementDirectionalXMonotoneTraits_2::ConstructOpposite_2` 

*/

class ArrDirectionalTraits::ConstructOpposite_2 {
public:

/// \name Has Models 
/// @{

/*! 
accepts an \f$ x\f$-monotone curve `xc` and returns its opposite curve, 
namely a curve whose graph is the same as `xc`'s, and whose source and 
target are swapped with respect to `xc`'s source and target. 
*/ 
ArrDirectionalTraits::X_monotone_curve_2 operator()(ArrDirectionalTraits::X_monotone_curve_2 xc); 

/// @}

}; /* end ArrDirectionalTraits::ConstructOpposite_2 */


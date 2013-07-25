namespace ArrDirectionalTraits {
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines `AdaptableBinaryFunction` 

\cgalHasModel `ArrangementDirectionalXMonotoneTraits_2::Intersect_2` 

*/

class Intersect_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
computes the intersections of `xc1` and `xc2` and 
inserts them <I>in an ascending lexicographic \f$ xy\f$-order</I> into the 
output iterator `oi`. The value-type of `Output_iterator` is 
`CGAL::Object`, where each `Object` wraps either a 
`pair<ArrDirectionalTraits::Point_2, ArrDirectionalTraits::Multiplicity>` object, which 
represents an intersection point with its multiplicity (in case the 
multiplicity is undefined or unknown, it is set to \f$ 0\f$) or an 
`ArrDirectionalTraits::X_monotone_curve_2` object, representing an 
overlapping subcurve of `xc1` and `xc2`. In the latter case, 
the overlapping subcurves are given the direction of `xc1` and 
`xc2` if their directions are identical. Otherwise, the overlapping 
subcurves are given an arbitrary direction. The operator returns a 
past-the-end iterator for the output sequence. 
*/ 
Output_iterator operator()(ArrDirectionalTraits::X_monotone_curve_2 xc1, 
ArrDirectionalTraits::X_monotone_curve_2 xc2, 
Output_iterator& oi); 

/// @}

}; /* end ArrDirectionalTraits::Intersect_2 */

}

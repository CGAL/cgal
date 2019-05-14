namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor 

\cgalHasModel ArrangementXMonotoneTraits_2::Intersect_2 

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
`pair<ArrTraits::Point_2,ArrTraits::Multiplicity>` object, which 
represents an intersection point with its multiplicity (in case the 
multiplicity is undefined or unknown, it should be set to \f$ 0\f$) or an 
`ArrTraits::X_monotone_curve_2` object, representing an 
overlapping subcurve of `xc1` and `xc2`. The operator 
returns a past-the-end iterator for the output sequence. 
*/ 
Output_iterator operator()(ArrTraits::X_monotone_curve_2 xc1, 
ArrTraits::X_monotone_curve_2 xc2, 
Output_iterator& oi); 

/// @}

}; /* end ArrTraits::Intersect_2 */

}

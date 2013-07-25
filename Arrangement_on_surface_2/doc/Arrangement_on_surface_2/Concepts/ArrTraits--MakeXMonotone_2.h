namespace ArrTraits {
/*!
\ingroup PkgArrangement2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor 

\cgalHasModel ArrangementTraits_2::Make_x_monotone_2 

*/

class MakeXMonotone_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
subdivides the input curve `c` into \f$ x\f$-monotone subcurves and 
isolated points, and inserts the results into a container through the 
given output iterator. The value type of `OutputIterator` is 
`CGAL::Object`, where each `Object` wraps either an 
`ArrTraits::X_monotone_curve_2` object or a `ArrTraits::Point_2` 
object. The operator returns a past-the-end iterator for the output 
sequence. 
*/ 
template <typename OutputIterator> 
OutputIterator operator()( ArrTraits::Curve_2 c, 
OutputIterator oi); 

/// @}

}; /* end ArrTraits::MakeXMonotone_2 */

}

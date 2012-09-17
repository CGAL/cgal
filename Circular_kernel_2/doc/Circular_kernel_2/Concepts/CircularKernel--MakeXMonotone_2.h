
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

For the sake of completeness, the `operator()` must also be 
defined for a 
`Line_arc_2`. In this case, the input line arc itself is the only 
arc returned through the `OutputIterator`. 

\sa `CircularKernel::MakeXYMonotone_2`

*/

class CircularKernel::MakeXMonotone_2 {
public:

/// \name See Also 
/// @{

/*! 
Splits the arc `ca` into \f$ x\f$-monotone arcs that are returned through the 
output iterator. Note that, to ensure an easy interface with the 
`Arrangement_2` package, the arcs are returned as `CGAL::Object`'s 
(see the `ArrangementTraits_2` concept). 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const CircularKernel::Circular_arc_2 &ca, OutputIterator oit); 

/// @}

}; /* end CircularKernel::MakeXMonotone_2 */



/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\sa `CircularKernel::MakeXYMonotone_2`

*/

class CircularKernel::MakeXMonotone_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Splits the arc `ca` into `x`-monotone arcs that are returned through the 
output iterator. Note that, to ensure an easy interface with the 
`Arrangement_2` package, the arcs are returned as `CGAL::Object`'s 
(see the `ArrangementTraits_2` concept). 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const CircularKernel::Circular_arc_2 &ca, OutputIterator oit); 

/// @}

}; /* end CircularKernel::MakeXMonotone_2 */


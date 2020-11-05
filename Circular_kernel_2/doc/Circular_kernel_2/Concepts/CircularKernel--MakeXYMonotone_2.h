
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\sa `CircularKernel::MakeXMonotone_2`

*/

class CircularKernel::MakeXYMonotone_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Splits the arc `ca` into `y`-monotone arcs that are returned through the
output iterator. Note that, to ensure an easy interface with the
`Arrangement_2` package, the arcs are returned as `CGAL::Object`'s
(see the `ArrangementTraits_2` concept).
*/
template < class OutputIterator >
OutputIterator
operator()(const CircularKernel::Circular_arc_2 &ca, OutputIterator oit);

/// @}

}; /* end CircularKernel::MakeXYMonotone_2 */


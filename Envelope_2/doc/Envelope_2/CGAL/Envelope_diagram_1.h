namespace CGAL {

/*!
\ingroup PkgEnvelope2

This class is the default envelope-diagram class used by envelope functions 
to represent the minimization or the maximization diagram of a set of curves. 
It represents the diagram as a doubly-linked list of interleaved vertices 
and edges. Thus, all operations provided by the envelope diagram take constant 
time, and the space needed to store the diagram class is linear in the 
complexity of the envelope. 

The envelope-diagram class is parameterized by a traits class, which is a 
model of the `ArrangementXMonotoneTraits_2` concept, in case we handle 
only envelopes of \f$ x\f$-monotone curves, or of the refined 
`ArrangementTraits_2` concept in case we handle arbitrary planar curves. 

\cgalModels `EnvelopeDiagram_1`

*/
template< typename Traits >
class Envelope_diagram_1 {
public:

/// @}

}; /* end Envelope_diagram_1 */
} /* end namespace CGAL */

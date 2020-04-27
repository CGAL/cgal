
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

Testing whether the interiors of two curves overlap.

*/

class CircularKernel::DoOverlap_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
For two line arcs.
*/
bool operator()
(const CircularKernel::Line_arc_2 & l0,
const CircularKernel::Line_arc_2 & l1);

/*!
For two circular arcs.
*/
bool operator()
(const CircularKernel::Circular_arc_2 & a0,
const CircularKernel::Circular_arc_2 & a1);

/// @}

}; /* end CircularKernel::DoOverlap_2 */



/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

To test whether a point lies in the vertical range of a curve.

*/

class CircularKernel::InXRange_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
For a line arc.
*/
bool operator()
(const CircularKernel::Line_arc_2 & l,
const CircularKernel::Circular_arc_point_2 & p);

/*!
For a circular arc. \pre `c` is `x`-monotone.
*/
bool operator()
(const CircularKernel::Circular_arc_2 & c,
const CircularKernel::Circular_arc_point_2 & p);

/// @}

}; /* end CircularKernel::InXRange_2 */


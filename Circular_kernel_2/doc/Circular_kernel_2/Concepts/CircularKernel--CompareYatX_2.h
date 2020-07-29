
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

A function object concept to compare a point `p` and an arc `a` on the vertical line passing through `p`.

*/

class CircularKernel::CompareYatX_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
For a circular arc.
\pre The arc `a` must be monotone and `p` must be in the vertical range of `a`.
*/
Comparison_result operator()
(const CircularKernel::Circular_arc_point_2 &p,
const CircularKernel::Circular_arc_2 &a);

/*!
Same for a segment.
*/
Comparison_result operator()
(const CircularKernel::Circular_arc_point_2 &p,
const CircularKernel::Line_arc_2 &a);

/// @}

}; /* end CircularKernel::CompareYatX_2 */


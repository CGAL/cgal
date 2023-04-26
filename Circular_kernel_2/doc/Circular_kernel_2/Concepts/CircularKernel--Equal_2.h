
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

Testing equality between objects.

\cgalRefines{Kernel::Equal_2}

\sa `CircularKernel::CompareX_2`
\sa `CircularKernel::CompareY_2`
\sa `CircularKernel::CompareXY_2`

*/

class CircularKernel::Equal_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
For two points.
*/
bool operator()
(const CircularKernel::Circular_arc_point_2 &p0,
const CircularKernel::Circular_arc_point_2 &p1);

/*!
For two arcs.
*/
bool operator()
(const CircularKernel::Circular_arc_2 &a0,
const CircularKernel::Circular_arc_2 &a1);

/*!
For two segments.
*/
bool operator()
(const CircularKernel::Line_arc_2 &a0,
const CircularKernel::Line_arc_2 &a1);

/// @}

}; /* end CircularKernel::Equal_2 */


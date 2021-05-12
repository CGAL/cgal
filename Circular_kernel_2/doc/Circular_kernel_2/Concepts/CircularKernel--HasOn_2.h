
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

To test whether a point lies on a curve.

\cgalRefines `Kernel::HasOn_2`
*/

class CircularKernel::HasOn_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
For a line.
*/
bool operator()
(const CircularKernel::Line_2 & l,
const CircularKernel::Circular_arc_point_2 &p);

/*!
For a circle.
*/
bool operator()
(const CircularKernel::Circle_2 & c,
const CircularKernel::Circular_arc_point_2 &p);

/*!
For a line arc.
*/
bool operator()
(const CircularKernel::Line_arc_2 & l,
const CircularKernel::Circular_arc_point_2 &p);

/*!
For a circular arc.
*/
bool operator()
(const CircularKernel::Circular_arc_2 & c,
const CircularKernel::Circular_arc_point_2 &p);

/// @}

}; /* end CircularKernel::HasOn_2 */


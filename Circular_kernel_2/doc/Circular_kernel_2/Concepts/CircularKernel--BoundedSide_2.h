
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::BoundedSide_2}

\sa `CircularKernel::HasOnBoundedSide_2`
\sa `CircularKernel::HasOnUnboundedSide_2`

*/

class CircularKernel::BoundedSide_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!

*/
Bounded_side
operator()
(const CircularKernel::Circle_2& s,
const CircularKernel::Circular_arc_point_2& p);

/// @}

}; /* end CircularKernel::BoundedSide_2 */


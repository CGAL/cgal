
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines `Kernel::CompareY_2`

\sa `CircularKernel::CompareX_2`
\sa `CircularKernel::CompareXY_2`
\sa `CircularKernel::Equal_2`

*/

class CircularKernel::CompareY_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the \f$ y\f$-coordinates of \f$ p\f$ and \f$ q\f$.
*/
Comparison_result operator()
(const CircularKernel::Circular_arc_point_2 &p,
const CircularKernel::Circular_arc_point_2 &q );

/// @}

}; /* end CircularKernel::CompareY_2 */


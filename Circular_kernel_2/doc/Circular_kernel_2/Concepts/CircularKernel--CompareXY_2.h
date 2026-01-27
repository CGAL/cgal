
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::CompareXY_2}

\sa `CircularKernel::CompareX_2`
\sa `CircularKernel::CompareY_2`
\sa `CircularKernel::Equal_2`

*/

class CircularKernel::CompareXY_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares \f$ p\f$ and \f$ q\f$ according to the lexicographic ordering on \f$ x\f$- and \f$ y\f$-coordinates.
*/
Comparison_result operator()
(const CircularKernel::Circular_arc_point_2 &p,
const CircularKernel::Circular_arc_point_2 &q );

/// @}

}; /* end CircularKernel::CompareXY_2 */


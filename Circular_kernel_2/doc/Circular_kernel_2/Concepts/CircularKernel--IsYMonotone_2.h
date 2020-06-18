
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::IsYMonotone_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Tests whether the arc is `y`-monotone.
*/
bool operator()
(const CircularKernel::Circular_arc_2 & c);

/*!
For a line arc, always returns `true`.
*/
bool operator()
(const CircularKernel::Line_arc_2 & l);

/// @}

}; /* end CircularKernel::IsYMonotone_2 */


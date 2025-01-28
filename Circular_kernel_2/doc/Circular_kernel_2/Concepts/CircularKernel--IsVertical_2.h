
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::IsVertical_2}
*/

class CircularKernel::IsVertical_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
For a line arc.
*/
bool operator()
(const CircularKernel::Line_arc_2 & l);

/*!
For a circular arc, always returns `false`.
*/
bool operator()
(const CircularKernel::Circular_arc_2 & c);

/// @}

}; /* end CircularKernel::IsVertical_2 */



/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::ConstructCircle_2}

\sa `CircularKernel::GetEquation`

*/

class CircularKernel::ConstructCircle_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs the supporting circle of a circular arc.
*/
CircularKernel::Circle_2 operator()
(CircularKernel::CircularArc_2);

/*!
Constructs a circle from an equation.
*/
CircularKernel::Circle_2 operator()
(CircularKernel::Polynomial_for_circles_2_2);

/// @}

}; /* end CircularKernel::ConstructCircle_2 */


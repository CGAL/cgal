
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\brief

\cgalRefines{Kernel::ConstructLine_3}

\sa `SphericalKernel::GetEquation`

*/

class SphericalKernel::ConstructLine_3 {
public:

/// \name Operations
///  A model of this concept must provide:
/// @{

/*!
Constructs the line containing the segment.
*/
SphericalKernel::Line_3 operator()
(const SphericalKernel::Line_arc_3 &s);

/*!
Constructs a line from an equation.
*/
SphericalKernel::Line_3 operator()
(const SphericalKernel::Polynomials_for_lines_3 &p);

/// @}

}; /* end SphericalKernel::ConstructLine_3 */


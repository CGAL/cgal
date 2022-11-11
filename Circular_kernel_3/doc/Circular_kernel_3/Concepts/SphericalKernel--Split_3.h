
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

*/

class SphericalKernel::Split_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Splits arc `a` at point `p`, which creates arcs `a1` and `a2`.
\pre The point `p` lies in the interior of the input arc `a`.
*/
void operator()
(const SphericalKernel::Circular_arc_3 &a,
const SphericalKernel::Circular_arc_point_3 &p,
SphericalKernel::Circular_arc_3 &a1,
SphericalKernel::Circular_arc_3 &a2);

/*!
Same for a line arc.
*/
void operator()
(const SphericalKernel::Line_arc_3 &l,
const SphericalKernel::Circular_arc_point_3 &p,
SphericalKernel::Line_arc_3 &l1, SphericalKernel::Line_arc_3 &l2);

/// @}

}; /* end SphericalKernel::Split_3 */


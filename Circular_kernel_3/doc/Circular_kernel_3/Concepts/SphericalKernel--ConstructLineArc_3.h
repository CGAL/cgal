
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

*/

class SphericalKernel::ConstructLineArc_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs the line segment supported by `l`, whose source
is `p` and whose target is `q`.
\pre `p` and `q` lie on `l` and are different.
*/
SphericalKernel::Line_arc_3 operator()
(const SphericalKernel::Line_3 &l,
const SphericalKernel::Circular_arc_point_3 &p,
const SphericalKernel::Circular_arc_point_3 &q);

/*!

*/
SphericalKernel::Line_arc_3 operator()
(const SphericalKernel::Segment_3 &s);

/*!

*/
SphericalKernel::Line_arc_3 operator()
(const SphericalKernel::Point_3 &p,
const SphericalKernel::Point_3 &q);

/// @}

}; /* end SphericalKernel::ConstructLineArc_3 */


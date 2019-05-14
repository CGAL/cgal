
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

Testing equality between objects.

\cgalRefines `Kernel::Equal_3`

\sa `SphericalKernel::CompareX_3`
\sa `SphericalKernel::CompareY_3`
\sa `SphericalKernel::CompareZ_3`
\sa `SphericalKernel::CompareXY_3`
\sa `SphericalKernel::CompareXYZ_3`

*/

class SphericalKernel::Equal_3 {
public:

/// \name Operations 
///  An object of this type must provide in addition:
/// @{

/*!
For two points. 
*/ 
bool operator() 
(const SphericalKernel::Circular_arc_point_3 &p0, 
const SphericalKernel::Circular_arc_point_3 &p1); 

/*!
For two arcs. Two arcs are equal, iff their non-oriented 
supporting planes are equal, and the centers and squared 
radii of their respective supporting circles are equal, and their 
sources and targets are equal. 
*/ 
bool operator() 
(const SphericalKernel::Circular_arc_3 &a0, 
const SphericalKernel::Circular_arc_3 &a1); 

/*!
For two segments. Two segments are equal, iff their non-oriented 
supporting lines are equal (i.e.\ they define the same set of 
points), and their endpoints are the same. 
*/ 
bool operator() 
(const SphericalKernel::Line_arc_3 &a0, 
const SphericalKernel::Line_arc_3 &a1); 

/// @}

}; /* end SphericalKernel::Equal_3 */


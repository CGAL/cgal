
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

*/

class SphericalKernel::ConstructCircularMaxVertex_3 {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!
Constructs the maximal vertex of `l` with lexicographically 
largest coordinates. 
*/ 
SphericalKernel::Circular_arc_point_3 operator() 
(const SphericalKernel::Line_arc_3 & l); 

/// @}

}; /* end SphericalKernel::ConstructCircularMaxVertex_3 */


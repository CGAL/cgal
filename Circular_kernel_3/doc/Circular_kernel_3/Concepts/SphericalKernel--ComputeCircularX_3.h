
/*!
\ingroup PkgSphericalKernel3Concepts
\cgalconcept

*/
class SphericalKernel::ComputeCircularX_3 {
public:
/*! 
returns the \f$ x\f$-coordinate of the point `p`. 
*/ 
SphericalKernel::Root_of_2 
operator()(const SphericalKernel::Circular_arc_point_3 & p) const; 

}; /* end SphericalKernel::ComputeCircularX_3 */



/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalconcept

*/

class SphericalKernel::Split_3 {
public:

/// \name Operations
/// A model `fo` of this type must provide: 
/// @{

/*! 
Splits arc \f$ a\f$ at point \f$ p\f$, which creates arcs \f$ a1\f$ and \f$ a2\f$. 
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


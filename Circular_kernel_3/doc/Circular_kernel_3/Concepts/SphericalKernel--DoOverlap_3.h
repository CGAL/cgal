
/*!
\ingroup PkgSphericalKernel3GeometricConcepts
\cgalConcept

*/
class SphericalKernel::DoOverlap_3 {
public:

/// \name Operations 
/// An object of this type must provide: 
/// @{

/*!
For two line arcs. The computation may be faster when the boolean is set to true.
*/ 
bool operator() 
(const SphericalKernel::Line_arc_3 & l0, 
const SphericalKernel::Line_arc_3 & l1, 
const bool known_equal_supporting_line = false); 

/*!
For two circular arcs. The computation may be faster when the boolean is set to true.
*/ 
bool operator() 
(const SphericalKernel::Circular_arc_3 & a0, 
const SphericalKernel::Circular_arc_3 & a1, 
const bool known_equal_supporting_circle = false); 

/// @}

}; /* end SphericalKernel::DoOverlap_3 */


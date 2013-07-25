
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Point_of_sphere_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the ith point defining the sphere \f$ s\f$. 
*/ 
bool operator()(const Kernel_d::Sphere_d& s, int i); 

/// @}

}; /* end Kernel_d::Point_of_sphere_d */


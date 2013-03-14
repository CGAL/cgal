
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Point_of_sphere_d {
public:

/*! 
returns the ith point defining the sphere \f$ s\f$. 
*/ 
bool operator()(const Kernel_d::Sphere_d& s, int i); 

}; /* end Kernel_d::Point_of_sphere_d */


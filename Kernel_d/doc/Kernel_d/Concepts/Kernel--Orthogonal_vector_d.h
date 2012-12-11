
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Orthogonal_vector_d {
public:

/*! 
computes an orthogonal vector to \f$ h\f$. 
*/ 
Kernel_d::Vector_d operator()(const Kernel_d::Hyperplane_d& h); 

}; /* end Kernel_d::Orthogonal_vector_d */


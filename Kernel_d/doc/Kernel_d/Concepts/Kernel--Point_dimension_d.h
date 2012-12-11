
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Point_dimension_d {
public:

/*! 
returns the dimension of \f$ p\f$. 
*/ 
int operator()(const Kernel_d::Point_d& 
p); 

}; /* end Kernel_d::Point_dimension_d */


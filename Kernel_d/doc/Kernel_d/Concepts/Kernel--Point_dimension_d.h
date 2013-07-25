
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Point_dimension_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the dimension of \f$ p\f$. 
*/ 
int operator()(const Kernel_d::Point_d& 
p); 

/// @}

}; /* end Kernel_d::Point_dimension_d */


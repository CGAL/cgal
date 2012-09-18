
/*!
\ingroup PkgKernelDKernelConcept
\cgalconcept

A model for this must provide: 

*/

class Kernel_d::Point_to_vector_d {
public:

/// \name See Also 
/// @{

/*! 
converts \f$ p\f$ to its geometric vector. 
*/ 
Kernel_d::Vector_d operator()(const Kernel_d::Point_d& p); 

/// @}

}; /* end Kernel_d::Point_to_vector_d */


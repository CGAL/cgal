
/*!
\ingroup PkgKernelDKernelConcept
\cgalconcept

A model for this must provide: 

*/

class Kernel_d::Value_at_d {
public:

/// \name See Also 
/// @{

/*! 
computes the value of \f$ h\f$ evaluated 
at \f$ p\f$.

\pre `p` and `h` have the same dimension. 
*/ 
Kernel_d::FT operator()(const Kernel_d::Hyperplane_d& 
h, const Kernel_d::Point_d& p); 

/// @}

}; /* end Kernel_d::Value_at_d */


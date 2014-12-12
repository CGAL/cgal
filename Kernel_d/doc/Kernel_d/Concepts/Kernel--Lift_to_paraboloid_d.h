
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class KernelWithLifting_d::Lift_to_paraboloid_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns \f$ p = (x_0,\ldots,x_{d-1})\f$ lifted to the paraboloid of 
revolution which is the point \f$ (x_0, \ldots,x_{d-1},\sum_{0 \le i <
d}x_i^2)\f$ in \f$ (d+1)\f$-space.
*/ 
Kernel_d::Point_d operator()(const Kernel_d::Point_d& 
p); 

/// @}

}; /* end KernelWithLifting_d::Lift_to_paraboloid_d */


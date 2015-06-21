
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class KernelWithLifting_d::Project_along_d_axis_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns \f$ p\f$ projected along the \f$ d\f$-axis onto the hyperspace 
spanned by the first \f$ d-1\f$ standard base vectors. 
*/ 
Kernel_d::Point_d operator()(const Kernel_d::Point_d& p); 

/// @}

}; /* end KernelWithLifting_d::Project_along_d_axis_d */


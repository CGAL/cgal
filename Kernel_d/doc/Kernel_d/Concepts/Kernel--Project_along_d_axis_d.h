
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Project_along_d_axis_d {
public:

/*! 
returns \f$ p\f$ projected along the \f$ d\f$-axis onto the hyperspace 
spanned by the first \f$ d-1\f$ standard base vectors. 
*/ 
Kernel_d::Point_d operator()(const Kernel_d::Point_d& p); 

}; /* end Kernel_d::Project_along_d_axis_d */


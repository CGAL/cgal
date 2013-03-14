
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Midpoint_d {
public:

/*! 
computes the midpoint of the segment 
\f$ pq\f$.\pre `p` and `q` have the same dimension. 
*/ 
Kernel_d::Point_d operator()(const Kernel_d::Point_d& 
p, const Kernel_d::Point_d& q); 

}; /* end Kernel_d::Midpoint_d */


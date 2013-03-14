
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Less_coordinate_d {
public:

/*! 
returns `true` iff the \f$ i\f$th %Cartesian coordinate 
of `p` is 
smaller than the \f$ i\f$th %Cartesian coordinate of `q`.

\pre `p` and `q` have the same dimension. 
*/ 
bool operator()(const Kernel_d::Point_d& 
p,const Kernel_d::Point_d& 
q, int i); 

}; /* end Kernel_d::Less_coordinate_d */



/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Less_coordinate_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns `true` iff the \f$ i\f$th %Cartesian coordinate 
of `p` is 
smaller than the \f$ i\f$th %Cartesian coordinate of `q`.

\pre `p` and `q` have the same dimension. 
*/ 
bool operator()(const Kernel_d::Point_d& 
p,const Kernel_d::Point_d& 
q, int i); 

/// @}

}; /* end Kernel_d::Less_coordinate_d */



/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Compare_lexicographically_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Compares the %Cartesian coordinates of 
points `p` and `q` lexicographically in ascending 
order of its %Cartesian components `p[i]` and `q[i]` for \f$ i = 
0,\ldots,d-1\f$.

\pre The objects are of the same dimension. 
*/ 
Comparison_result operator()(const Kernel_d::Point_d& 
p, const Kernel_d::Point_d& q); 

/// @}

}; /* end Kernel_d::Compare_lexicographically_d */



/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Orientation_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
determines the orientation of the points of the tuple 
`A = tuple [first,last)` where \f$ A\f$ consists of \f$ d + 1\f$ points in 
\f$ d\f$-space. This is the sign of the determinant 
\f[ \left| \begin{array}{cccc} 
1 & 1 & 1 & 1 \\ 
A[0] & A[1] & \dots& A[d] 
\end{array} \right| \f] 
where `A[i]` denotes the %Cartesian coordinate vector of 
the \f$ i\f$-th point in \f$ A\f$. 
\pre `size [first,last) == d+1` and `A[i].dimension() == d` \f$ \forall0 \leq i \leq d\f$. 

\tparam ForwardIterator has `Kernel_d::Point_d` as value type.
*/ 
template <class ForwardIterator> 
Orientation operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Orientation_d */


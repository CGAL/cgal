
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept
*/
class Kernel_d::Center_of_sphere_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the 
center of the sphere spanned by the points in `A = tuple [first,last)`.
\pre \f$A\f$ contains \f$d+1\f$ affinely independent points of dimension \f$d\f$. 
\tparam ForwardIterator has `Kernel_d::Point_d` as value type.
*/ 
template <class ForwardIterator> Kernel_d::Point_d 
operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Center_of_sphere_d */


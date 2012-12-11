
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Contained_in_affine_hull_d {
public:

/*! 
determines whether \f$ p\f$ is contained in the 
affine hull of the points in `A = tuple [first,last)`. 
\pre The objects are of the same dimension. 

\cgalRequires The value type of `ForwardIterator` is `Kernel_d::Point_d`. 
*/ 
template <class ForwardIterator> Bounded_side 
operator()( ForwardIterator first, ForwardIterator last, const 
Kernel_d::Point_d& p); 

}; /* end Kernel_d::Contained_in_affine_hull_d */


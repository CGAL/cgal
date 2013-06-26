
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Intersect_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! 
returns the result of the intersection of \f$ p\f$ and \f$ q\f$ in form of a 
polymorphic object. `Kernel_object` may be any of 
`Kernel_d::Segment_d`, `Kernel_d::Ray_d`, `Kernel_d::Line_d`, 
`Kernel_d::Hyperplane_d`. 

\pre `p` and `q` have the same dimension. 
*/ 
template <class Kernel_object> Object 
operator()(const Kernel_object& p, const Kernel_object& q); 

/// @}

}; /* end Kernel_d::Intersect_d */


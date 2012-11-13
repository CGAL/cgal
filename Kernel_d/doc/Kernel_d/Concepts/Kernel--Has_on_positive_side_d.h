
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Has_on_positive_side_d {
public:

/// \name See Also 
/// @{

/*! 
returns true iff \f$ p\f$ 
is on the positive side of \f$ o\f$. `Kernel_object` may be any of 
`Kernel_d::Sphere_d`, `Kernel_d::Hyperplane_d`.\pre `p` and `o` have the same dimension. 
*/ 
template <class Kernel_object> bool operator()(const 
Kernel_object& o, const Kernel_d::Point_d& p); 

/// @}

}; /* end Kernel_d::Has_on_positive_side_d */


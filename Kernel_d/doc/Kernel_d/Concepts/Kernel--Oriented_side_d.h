
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Oriented_side_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the side of \f$ p\f$ with respect to \f$ o\f$. `Kernel_object`
may be any of `Kernel_d::Sphere_d` or `Kernel_d::Hyperplane_d`.
\pre `p` and `o` have the same dimension.
*/
template <class Kernel_object> Oriented_side
operator()(const Kernel_object& o, const Kernel_d::Point_d& p);

/// @}

}; /* end Kernel_d::Oriented_side_d */


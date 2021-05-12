
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Value_at_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
computes the value of \f$ h\f$ evaluated
at \f$ p\f$.

\pre `p` and `h` have the same dimension.
*/
Kernel_d::FT operator()(const Kernel_d::Hyperplane_d&
h, const Kernel_d::Point_d& p);

/// @}

}; /* end Kernel_d::Value_at_d */


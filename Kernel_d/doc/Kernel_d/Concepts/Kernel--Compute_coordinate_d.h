
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Compute_coordinate_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the \f$ i\f$th cartesian coordinate of \f$ p\f$
*/
Kernel_d::FT operator()(const Kernel_d::Point_d&
p, int i);

/// @}

}; /* end Kernel_d::Compute_coordinate_d */



/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Vector_to_point_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
converts \f$ v\f$ to the affine point \f$ 0+v\f$.
*/
Kernel_d::Point_d operator()(const Kernel_d::Vector_d& v);

/// @}

}; /* end Kernel_d::Vector_to_point_d */


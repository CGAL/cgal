
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Point_to_vector_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
converts \f$ p\f$ to its geometric vector.
*/
Kernel_d::Vector_d operator()(const Kernel_d::Point_d& p);

/// @}

}; /* end Kernel_d::Point_to_vector_d */


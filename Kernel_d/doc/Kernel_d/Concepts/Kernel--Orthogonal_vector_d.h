
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Orthogonal_vector_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
computes an orthogonal vector to \f$ h\f$.
*/
Kernel_d::Vector_d operator()(const Kernel_d::Hyperplane_d& h);

/// @}

}; /* end Kernel_d::Orthogonal_vector_d */


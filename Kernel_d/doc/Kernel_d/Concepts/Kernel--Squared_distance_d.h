
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Squared_distance_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
computes the square of the Euclidean
distance between the two points \f$ p\f$ and \f$ q\f$.

\pre The dimensions of \f$ p\f$ and \f$ q\f$ are the same.
*/
Kernel_d::FT operator()(const Kernel_d::Point_d& p,
const Kernel_d::Point_d& q);

/// @}

}; /* end Kernel_d::Squared_distance_d */


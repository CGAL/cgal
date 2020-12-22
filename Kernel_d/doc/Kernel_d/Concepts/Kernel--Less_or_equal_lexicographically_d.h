
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Less_or_equal_lexicographically_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns `true` iff \f$ p\f$ is
lexicographically smaller than \f$ q\f$ with respect to %Cartesian
lexicographic order of points or equal to \f$ q\f$.

\pre `p` and `q` have the same dimension.
*/
bool operator()(const Kernel_d::Point_d& p, const
Kernel_d::Point_d& q);

/// @}

}; /* end Kernel_d::Less_or_equal_lexicographically_d */


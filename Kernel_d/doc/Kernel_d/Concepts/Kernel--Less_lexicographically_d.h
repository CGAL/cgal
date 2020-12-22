
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Less_lexicographically_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns `true` iff `p` is
lexicographically smaller than `q` with respect to %Cartesian
lexicographic order of points.

\pre `p` and `q` have the same dimension.
*/
bool operator()(const Kernel_d::Point_d&p, const
Kernel_d::Point_d&q);

/// @}

}; /* end Kernel_d::Less_lexicographically_d */


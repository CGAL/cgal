
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Side_of_oriented_sphere_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the relative position of point
`p` to the oriented sphere defined by the points in `A = tuple [first,last)` The order of the points in \f$ A\f$ is important,
since it determines the orientation of the implicitly constructed
sphere. If the points in \f$ A\f$ are positively oriented, the positive
side is the bounded interior of the sphere.

\pre `A` contains \f$ d+1\f$ points in \f$ d\f$-space.
\tparam ForwardIterator has `Kernel_d::Point_d` as value type.
*/
template <class ForwardIterator> Oriented_side
operator()( ForwardIterator first, ForwardIterator last, const
Kernel_d::Point_d& p);

/// @}

}; /* end Kernel_d::Side_of_oriented_sphere_d */


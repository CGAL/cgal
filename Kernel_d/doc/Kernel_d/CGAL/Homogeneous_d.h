
namespace CGAL {

/*!
\ingroup PkgKernelDKernels

A model for a `Kernel_d` (and even `KernelWithLifting_d`) using homogeneous coordinates to represent the
geometric objects. In order for `Homogeneous` to model Euclidean geometry
in \f$ E^d\f$, for some mathematical ring \f$ E\f$ (<I>e.g.</I>,
the integers \f$\mathbb{Z}\f$ or the rationals \f$\mathbb{Q}\f$), the template parameter `RT`
must model the mathematical ring \f$ E\f$. That is, the ring operations on this
number type must compute the mathematically correct results. If the number
type provided as a model for `RingNumberType` is only an approximation of a
ring (such as the built-in type `double`), then the geometry provided by
the kernel is only an approximation of Euclidean geometry.

\cgalModels `KernelWithLifting_d`

\sa `CGAL::Cartesian_d<FieldNumberType>`

*/
template< typename RingNumberType >
class Homogeneous_d {
public:

}; /* end Homogeneous */
} /* end namespace CGAL */

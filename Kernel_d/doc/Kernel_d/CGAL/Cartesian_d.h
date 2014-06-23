
namespace CGAL {

/*!
\ingroup PkgKernelDKernels

A model for `Kernel_d` (and even `KernelWithLifting_d`) that uses %Cartesian coordinates to represent the
geometric objects. In order for `Cartesian_d` to model Euclidean geometry 
in \f$ E^d\f$ , for some mathematical field \f$ E\f$ (<I>e.g.</I>, 
the rationals \f$\mathbb{Q}\f$ or the reals \f$\mathbb{R}\f$), the template parameter `FieldNumberType` 
must model the mathematical field \f$ E\f$. That is, the field operations on this 
number type must compute the mathematically correct results. If the number 
type provided as a model for `FieldNumberType` is only an approximation of a 
field (such as the built-in type `double`), then the geometry provided by 
the kernel is only an approximation of Euclidean geometry. 

\cgalModels `KernelWithLifting_d`

\sa `CGAL::Homogeneous_d<RingNumberType>`

*/
template< typename FieldNumberType >
class Cartesian_d {
public:

/// @}

}; /* end Cartesian_d */
} /* end namespace CGAL */

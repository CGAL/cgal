
namespace CGAL {

/*!
\ingroup kernel_classes

A model for a `Kernel` using homogeneous coordinates to represent the
geometric objects. In order for `Homogeneous` to model Euclidean geometry
in \f$ E^2\f$ and/or \f$ E^3\f$, for some mathematical ring \f$ E\f$ (<I>e.g.</I>,
the integers \f$\mathbb{Z}\f$ or the rationals \f$\mathbb{Q}\f$), the template parameter `RingNumberType`
must model the mathematical ring \f$ E\f$. That is, the ring operations on this
number type must compute the mathematically correct results. If the number
type provided as a model for `RingNumberType` is only an approximation of a
ring (such as the built-in type `double`), then the geometry provided by
the kernel is only an approximation of Euclidean geometry.

\cgalModels `Kernel`

\cgalHeading{Implementation}

This model of a kernel uses reference counting.

\sa `CGAL::Cartesian<FieldNumberType>`
\sa `CGAL::Simple_cartesian<FieldNumberType>`
\sa `CGAL::Simple_homogeneous<RingNumberType>`

*/
template< typename RingNumberType >
struct Homogeneous {

/// \name Types
/// @{

/*!

*/
typedef Quotient<RingNumberType> FT;

/*!

*/
typedef RingNumberType RT;

/// @}

}; /* end Homogeneous */
} /* end namespace CGAL */

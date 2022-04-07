
namespace CGAL {

/*!
\ingroup kernel_classes

A model for a `Kernel` using homogeneous coordinates to represent the
geometric objects. In order for `Simple_homogeneous` to model Euclidean geometry
in \f$ E^2\f$ and/or \f$ E^3\f$, for some mathematical ring \f$ E\f$ (<I>e.g.</I>,
the integers \f$\mathbb{Z}\f$ or the rationals \f$\mathbb{Q}\f$), the template parameter `RingNumberType`
must model the mathematical ring \f$ E\f$. That is, the ring operations on this
number type must compute the mathematically correct results. If the number
type provided as a model for `RingNumberType` is only an approximation of a
ring (such as the built-in type `double`), then the geometry provided by
the kernel is only an approximation of Euclidean geometry.

\cgalModels `Kernel`

\cgalHeading{Implementation}

In contrast to `Homogeneous`, no reference counting
is used internally. This eases debugging, but may slow down algorithms
that copy objects intensively, or slightly speed up others.

\sa `CGAL::Cartesian<FieldNumberType>`
\sa `CGAL::Homogeneous<RingNumberType>`
\sa `CGAL::Simple_cartesian<FieldNumberType>`

*/
template< typename RingNumberType >
struct Simple_homogeneous {

/// \name Types
/// @{

/*!

*/
typedef Quotient<RingNumberType> FT;

/*!

*/
typedef RingNumberType RT;

/// @}

}; /* end Simple_homogeneous */
} /* end namespace CGAL */

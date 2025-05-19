
namespace CGAL {

/*!
\ingroup kernel_classes

A model for `Kernel` that uses %Cartesian coordinates to represent the
geometric objects. In order for `Cartesian` to model Euclidean geometry
in \f$ E^2\f$ and/or \f$ E^3\f$, for some mathematical field \f$ E\f$ (<I>e.g.</I>,
the rationals \f$\mathbb{Q}\f$ or the reals \f$\mathbb{R}\f$), the template parameter `FieldNumberType`
must model the mathematical field \f$ E\f$. That is, the field operations on this
number type must compute the mathematically correct results. If the number
type provided as a model for `FieldNumberType` is only an approximation of a
field (such as the built-in type `double`), then the geometry provided by
the kernel is only an approximation of Euclidean geometry.

\cgalModels{Kernel}

\cgalHeading{Implementation}

All geometric objects in `Cartesian` are reference counted.

\sa `CGAL::Simple_cartesian<FieldNumberType>`
\sa `CGAL::Homogeneous<RingNumberType>`
\sa `CGAL::Simple_homogeneous<RingNumberType>`

*/
template< typename FieldNumberType >
struct Cartesian {

/// \name Types
/// @{

/*!

*/
typedef FieldNumberType FT;

/*!

*/
typedef FieldNumberType RT;

/// @}

}; /* end Cartesian */
} /* end namespace CGAL */

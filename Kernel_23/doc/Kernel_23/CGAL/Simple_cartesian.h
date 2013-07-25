
namespace CGAL {

/*!
\ingroup kernel_classes

A model for a `Kernel` using %Cartesian coordinates to represent the 
geometric objects. In order for `Simple_cartesian` to model Euclidean geometry 
in \f$ E^2\f$ and/or \f$ E^3\f$, for some mathematical field \f$ E\f$ (<I>e.g.</I>, 
the rationals \f$\mathbb{Q}\f$ or the reals \f$\mathbb{R}\f$), the template parameter `FieldNumberType` 
must model the mathematical field \f$ E\f$. That is, the field operations on this 
number type must compute the mathematically correct results. If the number 
type provided as a model for `FieldNumberType` is only an approximation of a 
field (such as the built-in type `double`), then the geometry provided by 
the kernel is only an approximation of Euclidean geometry. 

\cgalModels `Kernel`

\cgalHeading{Implementation}

In contrast to `Cartesian`, no reference counting 
is used internally. This eases debugging, but may slow down algorithms 
that copy objects intensively. 

\sa `CGAL::Cartesian<FieldNumberType>` 
\sa `CGAL::Homogeneous<RingNumberType>` 
\sa `CGAL::Simple_homogeneous<RingNumberType>` 

*/
template< typename FieldNumberType >
class Simple_cartesian {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef FieldNumberType FT; 

/*!

*/ 
typedef FieldNumberType RT; 

/// @}

}; /* end Simple_cartesian */
} /* end namespace CGAL */

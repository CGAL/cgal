
namespace CGAL {

/*!
\ingroup kernel_classes

The class `Kernel_traits` provides access to the kernel model to
which the argument type `T` belongs. (Provided `T` belongs to
some kernel model.) The default implementation assumes there is a
local type `T::R` referring to the kernel model of `T`.
If this type does not exist, a specialization of `Kernel_traits` can be
used to provide the desired information.

This class is, for example, useful in the following context. Assume
you want to write a generic function that accepts two points `p` and
`q` as argument and constructs the line segment between `p` and `q`.
In order to specify the return type of this function, you need to know
what is the segment type corresponding to the Point type representing
`p` and `q`. Using `Kernel_traits`, this can be done as follows.

\code
template < class Point >
typename Kernel_traits<Point>::Kernel::Segment
construct_segment(Point p, Point q)
{ ... }
\endcode

*/
template< typename T >
struct Kernel_traits {

/// \name Types
/// @{

/*!
If `T` is a type
`K::Point_2` of some kernel model `K`, then `Kernel` is
equal to `K`.
*/
typedef T::R Kernel;

/// @}

}; /* end Kernel_traits */
} /* end namespace CGAL */

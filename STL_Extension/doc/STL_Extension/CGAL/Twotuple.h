
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

The `Twotuple` class stores a homogeneous (same type) pair
of objects of type `T`. A `Twotuple` is much like a container, in that
it "owns" its elements. It is not actually a model of container, though,
because it does not support the standard methods (such as iterators) for
accessing the elements of a container.

\deprecated This class is deprecated, and will be removed in some future \cgal release.
Please use std::array instead.

\tparam T must be `Assignable`.

*/
template< typename T >
class Twotuple {
public:

/// \name Types
/// @{
/*!

*/
typedef T value_type;



/// @}


/// \name Variables
/// @{
/*!
first element
*/
T e0;



/// @}


/// \name Variables
/// @{
/*!
second element
*/
T e1;



/// @}


/// \name Creation
/// @{
/*!
introduces a `Twotuple` using the default
constructor of the elements.
*/
Twotuple();



/// @}


/// \name Creation
/// @{
/*!
constructs a `Twotuple` such
that `e0` is constructed from `x` and `e1` is
constructed from `y`.
*/
Twotuple(T x, T y);



/// @}



}; /* end Twotuple */
} /* end namespace CGAL */

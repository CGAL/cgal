
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities


\deprecated This class is deprecated, and will be removed in some future \cgal release.
Please use std::array instead.

The `Threetuple` class stores a homogeneous (same type) triple
of objects of type `T`. A `Threetuple` is much like a container, in that
it "owns" its elements. It is not actually a model of container, though,
because it does not support the standard methods (such as iterators) for
accessing the elements of a container.



\cgalHeading{Requirements}

`T` must be `Assignable`.


*/
template< typename T >
class Threetuple {
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


/// \name Variables
/// @{
/*!
third element
*/
T e2;



/// @}


/// \name Creation
/// @{
/*!
introduces a `Threetuple` using the default
constructor of the elements.
*/
Threetuple();



/// @}


/// \name Creation
/// @{
/*!
constructs a `Threetuple` such
that `e0` is constructed from `x`, `e1` is
constructed from `y` and `e2` is constructed from `z`.
*/
Threetuple(T x, T y, T z);



/// @}



}; /* end Threetuple */
} /* end namespace CGAL */

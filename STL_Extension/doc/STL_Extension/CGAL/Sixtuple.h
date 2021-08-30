
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

The `Sixtuple` class stores a homogeneous (same type)
sixtuple of objects of type `T`. A `Sixtuple` is much like a
container, in that it "owns" its elements. It is not actually a model of
container, though, because it does not support the standard methods (such as
iterators) for accessing the elements of a container.

\deprecated This class is deprecated, and will be removed in some future \cgal release.
Please use std::array instead.

\tparam T must be `Assignable`.
*/
template< typename T >
class Sixtuple {
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


/// \name Variables
/// @{
/*!
fourth element
*/
T e3;



/// @}


/// \name Variables
/// @{
/*!
fifth element
*/
T e4;



/// @}


/// \name Variables
/// @{
/*!
sixth element
*/
T e5;



/// @}


/// \name Creation
/// @{
/*!
introduces a `Sixtuple` using the default
constructor of the elements.
*/
Sixtuple();



/// @}


/// \name Creation
/// @{
/*!
constructs a
`Sixtuple` such that `e0` is constructed from `x`, `e1` from
`y`, `e2` from `z`, `e3` from `t`, `e4` from
`u` and `e5` from `v`.
*/
Sixtuple(T x, T y, T z, T t, T u, T v);



/// @}



}; /* end Sixtuple */
} /* end namespace CGAL */

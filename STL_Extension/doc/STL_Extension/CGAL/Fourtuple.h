
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities


The `Fourtuple` class stores a homogeneous (same type)
fourtuple of objects of type `T`. A `Fourtuple` is much like a
container, in that it "owns" its elements. It is not actually a model of
container, though, because it does not support the standard methods (such as
iterators) for accessing the elements of a container.

\deprecated This class is deprecated, and will be removed in some future \cgal release.
Please use std::array instead.

\tparam T must be `Assignable`.

*/
template< typename T >
class Fourtuple {
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


/// \name Creation
/// @{
/*!
introduces a `Fourtuple` using the default
constructor of the elements.
*/
Fourtuple();



/// @}


/// \name Creation
/// @{
/*!
constructs a `Fourtuple` such
that `e0` is constructed from `x`, `e1` from `y`,
`e2` from `z` and `e3` from `t`.
*/
Fourtuple(T x, T y, T z, T t);



/// @}



}; /* end Fourtuple */
} /* end namespace CGAL */

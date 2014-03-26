namespace CGAL {

/*!
\ingroup PkgPrincipalComponentAnalysisDbb

The function `bounding_box()` computes the axis-aligned bounding box of
a set of 2D or 3D points. The bounding box is returned either as an
iso rectangle in 2D or as an iso cuboid in 3D, the type being deduced
automatically from the value type of the iterator range.

There is a set of overloaded bounding_box functions for 2D and 3D
points. The user can also optionally pass an explicit kernel, in case
the default, based on `Kernel_traits` is not sufficient. The dimension
is also deduced automatically.
*/
/// @{

/*!
computes the bounding box of a non-empty set of 2D or 3D points.

\returns The return type is either `K::Iso_rectangle_2` or
`K::Iso_cuboid_3`, depending on the dimension of the input values,
where `K` is
\code
CGAL::Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel
\endcode

\pre first != beyond. 
*/
template < typename InputIterator >
Deduced
bounding_box(InputIterator first, InputIterator beyond);

/*!
computes the bounding box of a non-empty set of 2D or 3D points.

\returns The return type is either `K::Iso_rectangle_2` or
`K::Iso_cuboid_3`, depending on the dimension of the input values.

\pre first != beyond. 
*/
template < typename InputIterator, typename K >
Deduced
bounding_box(InputIterator first, InputIterator beyond, const K & k);

/// @}

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgPrincipalComponentAnalysisD

computes the bounding box of a non-empty set of 2D points.

The function `bounding_box` computes the axis-aligned bounding box of
a set of 2D or 3D points. The bounding box is returned either as an
iso rectangle in 2D or as an iso cuboid in 3D, the type being deduced
automatically from the value type of the iterator range.

The user can also optionally pass an explicit kernel, in case the
default, based on `Kernel_traits` is not sufficient. The dimension is
also deduced automatically.

`K` is `Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel`.
The value type must be `K::Point_2`.

\pre first != beyond. 
*/
template < typename InputIterator >
K::Iso_rectangle_2
bounding_box(InputIterator first, InputIterator beyond);

/*!
\ingroup PkgPrincipalComponentAnalysisD

computes the bounding box of a non-empty set of 2D points.

The function `bounding_box` computes the axis-aligned bounding box of
a set of 2D or 3D points. The bounding box is returned either as an
iso rectangle in 2D or as an iso cuboid in 3D, the type being deduced
automatically from the value type of the iterator range.

The user can also optionally pass an explicit kernel, in case the
default, based on `Kernel_traits` is not sufficient. The dimension is
also deduced automatically.


The value type must be `K::Point_2`.
\pre first != beyond. 
*/
template < typename InputIterator, typename K >
K::Iso_rectangle_2
bounding_box(InputIterator first, InputIterator beyond, const K & k);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the bounding box of a non-empty set of 3D points.

The function `bounding_box` computes the axis-aligned bounding box of
a set of 2D or 3D points. The bounding box is returned either as an
iso rectangle in 2D or as an iso cuboid in 3D, the type being deduced
automatically from the value type of the iterator range.

There is a set of overloaded `bounding_box` functions for 2D and 3D
points. The user can also optionally pass an explicit kernel, in case
the default, based on `Kernel_traits` is not sufficient. The dimension
is also deduced automatically.

`K` is `Kernel_traits<std::iterator_traits<InputIterator>::value_type>::Kernel`.
The value type must be `K::Point_3`.
\pre first != beyond. 
*/
template < typename InputIterator >
K::Iso_cuboid_3
bounding_box(InputIterator first, InputIterator beyond);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the bounding box of a non-empty set of 3D points.

The function `bounding_box` computes the axis-aligned bounding box of
a set of 2D or 3D points. The bounding box is returned either as an
iso rectangle in 2D or as an iso cuboid in 3D, the type being deduced
automatically from the value type of the iterator range.

There is a set of overloaded `bounding_box` functions for 2D and 3D
points. The user can also optionally pass an explicit kernel, in case
the default, based on `Kernel_traits` is not sufficient. The dimension
is also deduced automatically.

The value type must be `K::Point_3`.
\pre first != beyond. 
*/
template < typename InputIterator, typename K >
K::Iso_cuboid_3
bounding_box(InputIterator first, InputIterator beyond, const K & k);

} /* namespace CGAL */


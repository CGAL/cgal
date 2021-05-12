namespace CGAL {

/*!
\addtogroup PkgPrincipalComponentAnalysisDBary


The function `barycenter()` computes the barycenter (weighted center of
mass) of a set of 2D or 3D weighted points. The weight associated to
each point is specified using a `std::pair` storing the point and its
weight.

There is a set of overloaded `barycenter` functions for 2D and 3D
weighted points. The user can also optionally pass an explicit kernel,
in case the default, based on `Kernel_traits` is not sufficient. The
dimension is also deduced automatically.

\sa \link PkgPrincipalComponentAnalysisDCentroid `CGAL::centroid()` \endlink
\sa \link barycenter_grp `CGAL::barycenter() (Linear Kernel)` \endlink

*/
/// @{

/*!

computes the barycenter of a non-empty set of 2D or 3D weighted
points.

\returns `K::Point_2` or `K::Point_3` depending on the dimension of
the input values, where `K` is
\code
CGAL::Kernel_traits<
  std::iterator_traits<InputIterator>::value_type::first_type
>::Kernel
\endcode

\tparam InputIterator must have
`std::pair<K::Point_2, K::FT>` or `std::pair<K::Point_3, K::FT> as value type`.

\pre first != beyond, and the sum of the weights is non-zero.

*/
template < typename InputIterator >
Deduced
barycenter(InputIterator first, InputIterator beyond);

/*!
computes the barycenter of a non-empty set of 2D or 3D weighted
points.

\returns `K::Point_2` or `K::Point_3` depending on the dimension of
the input values.

\tparam InputIterator must have
`std::pair<K::Point_2, K::FT>` or `std::pair<K::Point_3, K::FT>` as value type.

\pre first != beyond, and the sum of the weights is non-zero.

*/
template < typename InputIterator, typename K >
Deduced
barycenter(InputIterator first, InputIterator beyond, const K & k);

/// @}
} /* namespace CGAL */


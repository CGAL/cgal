namespace CGAL {

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the barycenter of a non-empty set of 2D weighted points.

The function `barycenter` computes the barycenter (weighted center of
mass) of a set of 2D or 3D weighted points. The weight associated to
each point is specified using a `std::pair` storing the point and its
weight.

The user can also optionally pass an explicit kernel, in case the
default, based on `Kernel_traits` is not sufficient.  The dimension is
also deduced automatically.

\sa `CGAL::centroid` 


`K` is
`Kernel_traits<std::iterator_traits<InputIterator>::value_type::first_type>::Kernel`.
The value type must be `std::pair<K::Point_2, K::FT>`.

\pre first != beyond, and the sum of the weights is non-zero. 
*/
template < typename InputIterator >
K::Point_2
barycenter(InputIterator first, InputIterator beyond);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the barycenter of a non-empty set of 2D weighted points.

The function `barycenter` computes the barycenter (weighted center of
mass) of a set of 2D or 3D weighted points. The weight associated to
each point is specified using a `std::pair` storing the point and its
weight.

The user can also optionally pass an explicit kernel, in case the
default, based on `Kernel_traits` is not sufficient.  The dimension is
also deduced automatically.

\sa `CGAL::centroid` 

The value type must be `std::pair<K::Point_2, K::FT>`.

\pre first != beyond, and the sum of the weights is non-zero. 
*/
template < typename InputIterator, typename K >
K::Point_2
barycenter(InputIterator first, InputIterator beyond, const K & k);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the barycenter of a non-empty set of 3D weighted points.

The function `barycenter` computes the barycenter (weighted center of
mass) of a set of 2D or 3D weighted points. The weight associated to
each point is specified using a `std::pair` storing the point and its
weight.

The user can also optionally pass an explicit kernel, in case the
default, based on `Kernel_traits` is not sufficient.  The dimension is
also deduced automatically.

\sa `CGAL::centroid` 


`K` is
`Kernel_traits<std::iterator_traits<InputIterator>::value_type::first_type>::Kernel`.
The value type must be `std::pair<K::Point_3, K::FT>`.

\pre first != beyond, and the sum of the weights is non-zero. 
*/
template < typename InputIterator >
K::Point_3
barycenter(InputIterator first, InputIterator beyond);

/*!
\ingroup PkgPrincipalComponentAnalysisD

\brief computes the barycenter of a non-empty set of 3D weighted points.

The function `barycenter` computes the barycenter (weighted center of
mass) of a set of 2D or 3D weighted points. The weight associated to
each point is specified using a `std::pair` storing the point and its
weight.

The user can also optionally pass an explicit kernel, in case the
default, based on `Kernel_traits` is not sufficient.  The dimension is
also deduced automatically.

\sa `CGAL::centroid` 

The value type must be `std::pair<K::Point_3, K::FT>`.

\pre first != beyond, and the sum of the weights is non-zero. 
*/
template < typename InputIterator, typename K >
K::Point_3
barycenter(InputIterator first, InputIterator beyond, const K & k);

} /* namespace CGAL */


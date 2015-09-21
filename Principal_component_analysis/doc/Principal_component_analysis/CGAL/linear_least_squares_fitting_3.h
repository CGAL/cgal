namespace CGAL {

/*!
\ingroup PkgPrincipalComponentAnalysisDLLSF3

The function `linear_least_squares_fitting_3()` computes the best
fitting 3D line or plane (in the least squares sense) of a set of 3D
objects such as points, segments, triangles, spheres, balls, iso cuboids
or tetrahedra.

The best fitting linear sub-space (here line or plane) minimizes the
sum of squared distances from all points comprising these objects to
their orthogonal projections onto this linear subspace. It can be
shown that the best line or plane goes through the centroid of the
set. This problem is equivalent to search for the linear sub-space
which maximizes the variance of projected points (sum of squared
distances to the centroid). Internally we solve this problem by eigen
decomposition of the covariance matrix of the whole set. Note that the
\f$ 3 \times 3\f$ covariance matrix is computed internally in closed
form and not by point sampling the objects. Eigenvectors corresponding
to large eigenvalues are the directions in which the data has strong
component, or equivalently large variance.

The fitting quality property is characterized by the values of the
three eigenvalues. When all three values are distinct the best linear
subspace is uniquely determined, be it a line or a plane. When all
three eigenvalues are equal there is no preferable sub-space and any
line or plane going through the centroid share the same fitting
property (a horizontal plane or a line along the x axis are returned
by default). A best fitting line is uniquely determined as soon as the
largest eigenvalue is different from the two others, otherwise all
lines contained in the best fitting plane share the same fitting
property. A best fitting plane is uniquely determined as soon as the
smallest eigenvalue is different from the two others, otherwise all
planes going through the best fitting line share the same fitting
property.

*/
/// @{

/*!
\brief computes the best fitting 3D line of a 3D object set in the
range [`first`,`beyond`). The value returned is a fitting quality
between \f$ 0\f$ and \f$ 1\f$, where \f$ 0\f$ means that the variance
is the same along any line contained within the best fitting plane,
and \f$ 1\f$ means that the variance is null orthogonally to the best
fitting line (hence the fit is perfect).

The tag `tag` identifies the dimension to be considered from the
objects. For point sets it should be 0. For segment sets it could be 1
or 0 according to whether one wants to fit the entire segments or just
the end points. For triangle sets it can range from 0 to 2 according
to whether one wants to fit either the corner points, the segments or
the whole triangles. For iso cuboid sets it can range from 0 to 3
according to whether one wants to fit either the corners, the
segments, the faces or the whole solid iso cuboids. For sphere sets it can
be 2 or 3 according to whether one wants to fit either the surface of
the spheres or the whole solid balls. For tetrahedron sets it can
range from 0 to 3 according to whether one wants to fit either the
points, the segments, the surface triangles or the whole solid
tetrahedra.

The class `K` is the kernel in which the
value type of `InputIterator` is defined. It can be omitted and deduced
automatically from the value type.

The class `DiagonalizeTraits_` is a model of `DiagonalizeTraits`. It
can be omitted: if Eigen 3 (or greater) is available and
`CGAL_EIGEN3_ENABLED` is defined then an overload using
`Eigen_diagonalize_traits` is provided. Otherwise, the internal
implementation `Diagonalize_traits` is used.


\cgalHeading{Requirements}

<OL> 
<LI>`InputIterator` must have a value type  equivalent to `K::Point_3`, 
`K::Segment_3`, `K::Triangle_3`, `K::Iso_cuboid_3`, 
`K::Sphere_3` or `K::Tetrahedron_3`. 
<LI>`line` is the best fitting line computed. 
<LI>`centroid` is the centroid computed. This parameter is optional and can be omitted. 
<LI>`tag` is the tag identifying the dimension to be considered from the objects. It should range from `Dimension_tag<0>` to `Dimension_tag<3>`. Also, it should not be of a dimension greater nor smaller than the geometry of the object. For example, a `Triangle` can not have a `Dimension_tag<3>` tag. A `Segment` can not have a `Dimension_tag<2>` nor a `Dimension_tag<3>` tag. A `Sphere` can not have a `Dimension_tag<0>` nor a `Dimension_tag<1>` tag. 
</OL> 


*/
template < typename InputIterator, typename K, typename Tag, typename DiagonalizeTraits_ > 
typename K::FT 
linear_least_squares_fitting_3(InputIterator first,
InputIterator beyond,
typename K::Line_3& line,
typename K::Point_3& centroid, 
const Tag& tag,
const K& k,
const DiagonalizeTraits_& diagonalize_traits);

/*!
\brief computes the best fitting 3D plane of a 3D object set in the
range [`first`,`beyond`). The value returned is a fitting quality
between \f$ 0\f$ and \f$ 1\f$, where \f$ 0\f$ means that the variance
is the same along any plane going through the best fitting line, and
\f$ 1\f$ means that the variance is null orthogonally to the best
fitting plane (hence the fit is perfect).


The class `K` is the kernel in which the value type
of `InputIterator` is defined. It can be omitted and deduced
automatically from the value type. The tag `tag` identifies the
dimension to be considered from the objects (see above).

The class `DiagonalizeTraits_` is a model of `DiagonalizeTraits`. It
can be omitted: if Eigen 3 (or greater) is available and
`CGAL_EIGEN3_ENABLED` is defined then an overload using
`Eigen_diagonalize_traits` is provided. Otherwise, the internal
implementation `Diagonalize_traits` is used.

\cgalHeading{Requirements}

<OL> 
<LI>`InputIterator` has a value type equivalent to `K::Point_3`, 
`K::Segment_3`, `K::Triangle_3`, `K::Iso_cuboid_3`, 
`K::Sphere_3` or `K::Tetrahedron_3`. 
<LI>`plane` is the best fitting plane computed. 
<LI>`centroid` is the centroid computed. This parameter is optional and can be omitted. 
<LI>`tag` is the tag identifying the dimension to be considered from the objects. It should range from `Dimension_tag<0>` to `Dimension_tag<3>`. Also, it should not be of a dimension greater nor smaller than the geometry of the object. For example, a `Triangle` can not have a `Dimension_tag<3>` tag. A `Segment` can not have a `Dimension_tag<2>` nor a `Dimension_tag<3>` tag. A `Sphere` can not have a `Dimension_tag<0>` nor a `Dimension_tag<1>` tag. 
</OL> 

*/
template < typename InputIterator, typename K, typename Tag, typename DiagonalizeTraits_ > 
typename K::FT 
linear_least_squares_fitting_3(InputIterator first,
InputIterator beyond,
typename K::Plane_3& plane,
typename K::Point_3& centroid, 
const Tag& tag,
const K& k,
const DiagonalizeTraits_& diagonalize_traits);

/// @}

} /* namespace CGAL */


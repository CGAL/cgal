namespace CGAL {

/*!
\ingroup PkgPrincipalComponentAnalysisDLLSF2

\brief computes the best fitting 2D line of a 2D object set in the range [`first`,`beyond`). The value returned is a fitting quality between \f$ 0\f$ and \f$ 1\f$, where \f$ 0\f$ means that the variance is the same along any line (a horizontal line going through the centroid is output by default), and \f$ 1\f$ means that the variance is null orthogonally to the best fitting line (hence the fit is perfect).

It computes the 2D best fitting line (in the least squares sense) of a set of 2D objects such as points, segments, triangles, iso rectangles, circles or disks. 

The best fitting line minimizes the sum of squared distances from all points comprising these objects to their orthogonal projections onto the line. It can be shown that this line goes through the centroid of the set. This problem is equivalent to search for the linear sub-space which maximizes the variance of projected points (sum of squared distances to the centroid). Internally we solve this problem by eigen decomposition of the covariance matrix of the whole set. Note that the \f$ 2 \times 2\f$ covariance matrix is computed internally in closed form and not by point sampling the objects. Eigenvectors corresponding to large eigenvalues are the directions in which the data has strong component, or equivalently large variance. If one eigenvalue is null the fit is perfect as the sum of squared distance from all points to their projection onto the best line is null. If the two eigenvalues are the same there is no preferable sub-space and all lines going through the centroid share the same fitting property. 

The tag `tag` identifies the dimension to be considered from the objects. For point sets it should be 0. For segments it can be 1 or 0 according to whether one wants to fit the whole segment or just their end points. For triangles it can range from 0 to 2 according to whether one wants to fit either the triangle points, the segments or the whole triangles. For rectangles it can range from 0 to 2 according to whether one wants to fit either the corner points, the segments, or the whole rectangles. For circles it can be 1 or 2 according to whether one wants to fit either the circles or the whole discs. For triangles it ranges from 0 to 2 according to whether one wants to fit either the points, the segments or the whole triangles. 

The class `K` is the kernel in which the value type of the  `InputIterator` is defined. It can be omitted and deduced automatically from the value type. 

The class `DiagonalizeTraits_` is a model of `DiagonalizeTraits`. It can be omitted: if Eigen 3 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined then an overload using `Eigen_diagonalize_traits` is provided. Otherwise, the internal implementation `Diagonalize_traits` is used.

\cgalHeading{Requirements}

<OL> 
<LI>`InputIterator` must have a value type equivalent to `K::Point_2` or 
`K::Segment_2` or `K::Triangle_2` or `K::Rectangle_2` or 
`K::Circle_2`. 
<LI>`line` is the best fitting line computed. 
<LI>`centroid` is the centroid computed. This parameter is optional and can be 
omitted. 
<LI>`tag` is the tag identifying the dimension to be considered from the objects. It should be one of `Dimension_tag<0>`, `Dimension_tag<1>` or `Dimension_tag<2>`. Also, it should not be of dimension greater than the geometry of the object. For example, a `Segment` can not have a `Dimension_tag<2>` tag. 
</OL> 


\pre first != beyond. 
*/
template < typename InputIterator, typename K, typename Tag, typename DiagonalizeTraits_ > 
typename K::FT 
linear_least_squares_fitting_2(InputIterator first,
InputIterator beyond,
typename K::Line_2 & line,
typename K::Point_2 & centroid, 
const Tag & tag,
const K & k,
const DiagonalizeTraits_& diagonalize_traits);

} /* namespace CGAL */


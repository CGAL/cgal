namespace CGAL {

/*!
\ingroup DistanceClasses

The class `Weighted_Minkowski_distance` provides an implementation of the concept `OrthogonalDistance`, with a weighted
Minkowski metric on \f$ d\f$-dimensional points
defined by \f$ l_p(w)(r,q)= ({\Sigma_{i=1}^{i=d} \, w_i(r_i-q_i)^p})^{1/p}\f$ for \f$ 0 < p <\infty\f$ and
defined by \f$ l_{\infty}(w)(r,q)=max \{w_i |r_i-q_i| \mid 1 \leq i \leq d\}\f$.
For the purpose of the distance computations it is more efficient to compute
the transformed distance \f$ {\sigma_{i=1}^{i=d} \, w_i(r_i-q_i)^p}\f$ instead of the actual distance.

\tparam Traits must be a model of the concept
`SearchTraits`, for example `Search_traits_2`.

\cgalModels `OrthogonalDistance`

\sa `OrthogonalDistance`
\sa `CGAL::Euclidean_distance<Traits>`

*/
template< typename Traits >
class Weighted_Minkowski_distance {
public:

/// \name Types
/// @{

/*!
Dimension tag.
*/
typedef Traits::Dimension D;

/*!
Number type.
*/
typedef Traits::FT FT;

/*!
Point type.
*/
typedef Traits::Point_d Point_d;

/// @}

/// \name Creation
/// @{

/*!
Constructor implementing \f$ l_2\f$ metric for \f$ d\f$-dimensional points.
*/
Weighted_Minkowski_distance(int d,Traits t=Traits());

/*!
Constructor implementing the \f$ l_{power}(weights)\f$ metric. \f$ power \leq0\f$ denotes the \f$ l_{\infty}(weights)\f$ metric.
The values in the iterator range `[wb,we)` are the weight.
*/
template <class InputIterator>
Weighted_Minkowski_distance(FT power, int dim, InputIterator wb, InputIterator we,Traits t=Traits());

/// @}

/// \name Operations
/// @{

/*!
Returns \f$ d^{power}\f$,
where \f$ d\f$ denotes the distance between `q` and `r`.
*/
FT transformed_distance(Point_d q, Point_d r) const;

/*!
Returns \f$ d^{power}\f$, where \f$ d\f$ denotes the distance between the query item `q` and
the point on the boundary of `r` closest to `q`.
*/
FT min_distance_to_rectangle(Point_d q, Kd_tree_rectangle<FT,D> r) const;

/*!
Returns \f$ d^{power}\f$, where \f$ d\f$ denotes the distance between the query item `q` and
the point on the boundary of `r` closest to `q`. Stores the distances in each
dimension in `dists`.
*/
FT min_distance_to_rectangle(Point_d q, Kd_tree_rectangle<FT,D> r, vector<FT>& dists);

/*!
Returns \f$ d^{power}\f$, where \f$ d\f$ denotes the distance between the query item `q` and
the point on the boundary of `r` farthest to `q`.
*/
FT max_distance_to_rectangle(Point_d q, Kd_tree_rectangle<FT,D> r) const;

/*!
Returns \f$ d^{power}\f$, where \f$ d\f$ denotes the distance between the query item `q` and
the point on the boundary of `r` farthest to `q`. Stores the distances in each
dimension in `dists`.
*/
FT max_distance_to_rectangle(Point_d q, Kd_tree_rectangle<FT,D> r, vector<FT>& dists) ;

/*!
Updates `dist` incrementally
and returns the updated distance.
*/
FT new_distance(FT dist, FT old_off, FT new_off, int cutting_dimension) const;

/*!
Returns \f$ d^p\f$ for \f$ 0 < p <\infty\f$ . Returns \f$ d\f$ for \f$ p=\infty\f$ .
*/
FT transformed_distance(FT d) const;

/*!
Returns \f$ d^{1/p}\f$ for \f$ 0 < p <\infty\f$.
Returns \f$ d\f$ for \f$ p=\infty\f$.
*/
FT inverse_of_transformed_distance(FT d) const;

/// @}

}; /* end Weighted_Minkowski_distance */
} /* end namespace CGAL */

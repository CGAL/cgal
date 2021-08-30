namespace CGAL {

/// \addtogroup PkgKernelDFunctions
/// @{

/*!


returns the center of the sphere spanned by the points in `A = tuple[first,last)`.

\pre \f$ A\f$ contains \f$ d+1\f$ affinely independent points of dimension \f$ d\f$.
\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> Point_d<R>
center_of_sphere(ForwardIterator first, ForwardIterator last);


/*!


returns the projection of \f$ p = (x_0,\ldots,x_{d-1})\f$ onto the
paraboloid of revolution which is the point \f$ (p_0,
\ldots,p_{d-1},\sum_{0 \le i < d}p_i^2)\f$ in \f$ (d+1)\f$-space.
*/
Point_d<R> lift_to_paraboloid(const Point_d<R>& p);


/*!


computes a basis of the linear space
spanned by the vectors in `A = tuple [first,last)` and returns
it via an iterator range starting in `result`. The returned
iterator marks the end of the output.
\pre \f$ A\f$ contains vectors of the same dimension \f$ d\f$.

\tparam ForwardIterator has `Vector_d<R>` as value type.
*/
template <class ForwardIterator, class OutputIterator>
OutputIterator linear_base(ForwardIterator first, ForwardIterator
last, OutputIterator result);


/*!


computes the midpoint of the segment \f$ pq\f$.

\pre `p.dimension() == q.dimension()`.
*/
Point_d<R> midpoint(const Point_d<R>& p, const Point_d<R>&
q);


/*!


returns \f$ p\f$ projected along the \f$ d\f$-axis onto the hyperspace spanned
by the first \f$ d-1\f$ standard base vectors.
*/
Point_d<R> project_along_d_axis(const Point_d<R>& p);


/*!


computes the square of the Euclidean distance between the two points
\f$ p\f$ and \f$ q\f$.

\pre The dimensions of \f$ p\f$ and \f$ q\f$ are the same.
*/
FT squared_distance(Point_d<R> p, Point_d<R> q);

/// @}

} /* namespace CGAL */


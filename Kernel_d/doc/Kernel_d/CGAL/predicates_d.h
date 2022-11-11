namespace CGAL {

/// \addtogroup PkgKernelDFunctions
/// @{

/*!


returns true iff the points in `A = tuple [first,last)` are
affinely independent.

\pre The objects are of the same dimension.
\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> bool
affinely_independent(ForwardIterator first, ForwardIterator last);


/*!


computes
the affine rank of the points in `A = tuple [first,last)`.
\pre The objects in \f$ A\f$ are of the same dimension.
\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> int
affine_rank(ForwardIterator first, ForwardIterator last);


/*!


Compares the %Cartesian
coordinates of points `p` and `q` lexicographically
in ascending order of its %Cartesian components `p[i]` and
`q[i]` for \f$ i = 0,\ldots,d-1\f$.

\pre `p.dimension() == q.dimension()`.
*/
Comparison_result compare_lexicographically(const
Point_d<R>& p, const Point_d<R>& q);


/*!


determines whether \f$ p\f$ is contained in
the affine hull of the points in `A = tuple [first,last)`.
\pre The objects in `A` are of the same dimension.

\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> bool
contained_in_affine_hull( ForwardIterator first, ForwardIterator
last, const Point_d<R>& p);


/*!


determines whether \f$ v\f$ is contained
in the linear hull of the vectors in `A = tuple [first,last)`.
\pre The objects in \f$ A\f$ are of the same dimension.

\tparam ForwardIterator has `Vector_d<R>` as value type.
*/
template <class ForwardIterator> bool
contained_in_linear_hull( ForwardIterator first, ForwardIterator
last, const Vector_d<R>& v);


/*!


determines whether \f$ p\f$ is contained in the
simplex of the points in `A = tuple [first,last)`.

\pre The objects in \f$ A\f$ are of the same dimension and affinely
independent.
\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> bool
contained_in_simplex( ForwardIterator first, ForwardIterator last,
const Point_d<R>& p);


/*!


returns `true` iff `p` is
lexicographically smaller than `q` with respect to %Cartesian
lexicographic order of points.

\pre `p.dimension() == q.dimension()`.
*/
bool lexicographically_smaller(const Point_d<R>& p, const
Point_d<R>& q);


/*!


returns `true` iff \f$ p\f$ is
lexicographically smaller than \f$ q\f$ with respect to %Cartesian
lexicographic order of points or equal to \f$ q\f$.

\pre `p.dimension() == q.dimension()`.
*/
bool lexicographically_smaller_or_equal( const Point_d<R>&
p, const Point_d<R>& q);


/*!


decides whether the vectors in `A = tuple [first,last)`
are linearly independent.

\pre The objects in `A` are of the same dimension.
\tparam ForwardIterator has `Vector_d<R>` as value type.
*/
template <class ForwardIterator> bool
linearly_independent( ForwardIterator first, ForwardIterator
last);


/*!


computes
the linear rank of the vectors in `A = tuple [first,last)`.
\pre The objects are of the same dimension.
\tparam ForwardIterator has `Vector_d<R>` as value type.
*/
template <class ForwardIterator> int
linear_rank(ForwardIterator first, ForwardIterator last);


/*!


determines the orientation of the points of the tuple `A = tuple [first,last)` where \f$ A\f$ consists of \f$ d+1\f$ points in
\f$ d\f$-space. This is the sign of the determinant
\f[ \left| \begin{array}{cccc}
1 & 1 & 1 & 1 \\
A[0] & A[1] & \dots& A[d]
\end{array} \right| \f]
where `A[i]` denotes the %Cartesian coordinate vector of
the \f$ i\f$-th point in \f$ A\f$.
\pre `size [first,last) == d+1` and `A[i].dimension() == d` \f$ \forall0 \leq i \leq d\f$.

\tparam ForwardIterator has `Point_d<R>` as value type.

*/
template <class ForwardIterator> Orientation
orientation(ForwardIterator first, ForwardIterator last);


/*!


returns the relative position of point
`p` to the sphere defined by `A = tuple [first,last)`. The
order of the points of \f$ A\f$ does not matter.

\pre `orientation(first,last)` is not `ZERO`.
\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> Bounded_side
side_of_bounded_sphere( ForwardIterator first, ForwardIterator last,
const Point_d<R>& p);

/*!


returns the relative position of
point `p` to the oriented sphere defined by the points in
`A = tuple [first,last)` The order of the points in \f$ A\f$ is
important, since it determines the orientation of the implicitly
constructed sphere. If the points in \f$ A\f$ are positively oriented,
the positive side is the bounded interior of the sphere.

\pre `A` contains \f$ d+1\f$ points in \f$ d\f$-space.
\tparam ForwardIterator has `Point_d<R>` as value type.
*/
template <class ForwardIterator> Oriented_side
side_of_oriented_sphere( ForwardIterator first, ForwardIterator
last, const Point_d<R>& p);

/// @}

} /* namespace CGAL */


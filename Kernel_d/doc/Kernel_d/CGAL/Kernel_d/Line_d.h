namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance of data type `Line_d` is an oriented line in
\f$ d\f$-dimensional Euclidean space.

\cgalHeading{Implementation}

Lines are implemented by a pair of points as an item type. All
operations like creation, initialization, tests, direction
calculation, input and output on a line \f$ l\f$ take time
\cgalBigO{l.dimension()}. `dimension()`, coordinate and point
access, and identity test take constant time. The operations for
intersection calculation also take time \cgalBigO{l.dimension()}. The
space requirement is \cgalBigO{l.dimension()}.

*/
template< typename Kernel >
class Line_d {
public:

/// \name Types
/// @{

/*!
the linear algebra layer.

*/
typedef unspecified_type LA;

/// @}

/// \name Creation
/// @{

/*!
introduces a variable `l` of
type `Line_d<Kernel>`.
*/
Line_d<Kernel>();

/*!
introduces a
line through `p` and `q` and oriented from `p` to
`q`.

\pre \f$ p\f$ and \f$ q\f$ are distinct and have the same dimension.
*/
Line_d<Kernel>(Point_d<Kernel> p, Point_d<Kernel> q);

/*!
introduces
a line through `p` with direction `dir`.

\pre `p.dimension()==dir.dimension()`, `dir` is not degenerate.
*/
Line_d<Kernel>(Point_d<Kernel> p, Direction_d<Kernel> dir);

/*!
introduces a variable
`l` of type `Line_d<Kernel>` and initializes it to the line through
`s.source()` and `s.target()` with direction from
`s.source()` to `s.target()`.

\pre \f$ s\f$ is not degenerate.
*/
Line_d<Kernel>(Segment_d<Kernel> s);

/*!
introduces a variable `l` of
type `Line_d<Kernel>` and initializes it to the line through
`r.point(1)` and `r.point(2)`.
*/
Line_d<Kernel>(Ray_d<Kernel> r);

/// @}

/// \name Operations
/// @{

/*!
returns the dimension of the ambient
space.
*/
int dimension();

/*!
returns an arbitrary point on
`l`. It holds that `point(i) == point(j)`, iff
`i==j`. Furthermore, `l` is directed from `point(i)` to
`point(j)`, for all `i < j`.
*/
Point_d<Kernel> point(int i) ;

/*!
returns the line
`(point(2),point(1))` of opposite direction.
*/
Line_d<Kernel> opposite() ;

/*!
returns the direction of
`l`.
*/
Direction_d<Kernel> direction();

/*!
returns \f$ t(l)\f$.

\pre `l.dimension()==t.dimension()`.
*/
Line_d<Kernel> transform(const Aff_transformation_d<Kernel> &
t);

/*!
returns
`l+v`, i.e., `l` translated by vector \f$ v\f$.
\pre `l.dimension()==v.dimension()`.
*/
Line_d<Kernel> operator+(const Vector_d<Kernel>& v);

/*!
returns the
point of intersection of `l` with the hyperplane that is
orthogonal to `l` and that contains `p`.

\pre `l.dimension()==p.dimension()`.
*/
Point_d<Kernel> projection(const Point_d<Kernel>& p) ;

/*!
returns true if \f$ p\f$ lies
on `l` and false otherwise.

\pre `l.dimension()==p.dimension()`.
*/
bool has_on(const Point_d<Kernel>& p) ;

/// @}

}; /* end Line_d */


/*!
Test for equality as unoriented lines

\pre `l1.dimension()==l2.dimension()`.
\relates Line_d
*/
bool weak_equality(const Line_d<Kernel>& l1, const Line_d<Kernel>& l2) ;

/*!
returns true if `l1` and `l2` are parallel as unoriented
lines and false otherwise.

\pre `l1.dimension()==l2.dimension()`.
\relates Line_d
*/
bool parallel(const Line_d<Kernel>& l1, const Line_d<Kernel>& l2) ;
} /* end namespace CGAL */

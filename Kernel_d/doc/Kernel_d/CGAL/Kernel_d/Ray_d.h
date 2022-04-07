namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance of data type `Ray_d` is a ray in \f$ d\f$-dimensional
Euclidean space. It starts in a point called the source of `r` and
it goes to infinity.

\cgalHeading{Implementation}

Rays are implemented by a pair of points as an item type. All
operations like creation, initialization, tests, direction
calculation, input and output on a ray \f$ r\f$ take time
\f$ O(r.dimension())\f$. `dimension()`, coordinate and point
access, and identity test take constant time. The space requirement is
\f$ O(r.dimension())\f$.

*/
template< typename Kernel >
class Ray_d {
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
introduces some ray in
\f$ d\f$-dimensional space.
*/
Ray_d<Kernel>();

/*!
introduces a ray
through `p` and `q` and starting at `p`.

\pre \f$ p\f$ and \f$ q\f$ are distinct and have the same dimension.

\pre `p.dimension()==q.dimension()`.
*/
Ray_d<Kernel>(Point_d<Kernel> p, Point_d<Kernel> q);

/*!
introduces
a ray starting in `p` with direction `dir`.

\pre `p` and `dir` have the same dimension and `dir` is not degenerate.

\pre `p.dimension()==dir.dimension()`.
*/
Ray_d<Kernel>(Point_d<Kernel> p, Direction_d<Kernel> dir);

/*!
introduces a ray through
`s.source()` and `s.target()` and starting at
`s.source()`.

\pre \f$ s\f$ is not degenerate.
*/
Ray_d<Kernel>(Segment_d<Kernel> s);

/// @}

/// \name Operations
/// @{

/*!
returns the dimension of the ambient
space.
*/
int dimension() ;

/*!
returns the source point of `r`.

*/
Point_d<Kernel> source() ;

/*!
returns a point on
`r`. `point(0)` is the source. `point(i)`, with \f$ i>0\f$, is
different from the source.

\pre \f$ i \geq0\f$.
*/
Point_d<Kernel> point(int i) ;

/*!
returns the direction of
`r`.
*/
Direction_d<Kernel> direction() ;

/*!
returns the supporting line
of `r`.
*/
Line_d<Kernel> supporting_line() ;

/*!
returns the ray with direction
opposite to `r` and starting in `source`.
*/
Ray_d<Kernel> opposite() ;

/*!
returns \f$ t(r)\f$.

\pre `r.dimension()==t.dimension()`.
*/
Ray_d<Kernel> transform(const Aff_transformation_d<Kernel>& t)
;

/*!
returns
`r+v`, i.e., `r` translated by vector \f$ v\f$.

\pre `r.dimension()==v.dimension()`.
*/
Ray_d<Kernel> operator+(const Vector_d<Kernel>& v) ;

/*!
A point is on `r`,
iff it is equal to the source of `r`, or if it is in the interior
of `r`.

\pre `r.dimension()==p.dimension()`.
*/
bool has_on(const Point_d<Kernel>& p) ;


/// @}

}; /* end Ray_d */

/*!
returns true if the unoriented supporting lines of `r1` and
`r2` are parallel and false otherwise.

\pre `r1.dimension()==r2.dimension()`.
\relates Ray_d
*/
bool parallel(const Ray_d<Kernel>& r1, const Ray_d<Kernel>& r2) ;

} /* end namespace CGAL */

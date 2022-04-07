namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance \f$ S\f$ of the data type `Sphere_d` is an oriented sphere
in some \f$ d\f$-dimensional space. A sphere is defined by \f$ d+1\f$ points
(class `Point_d<Kernel>`). We use \f$ A\f$ to denote the array of the
defining points. A set \f$ A\f$ of defining points is <I>legal</I> if
either the points are affinely independent or if the points are all
equal. Only a legal set of points defines a sphere in the geometric
sense and hence many operations on spheres require the set of defining
points to be legal. The orientation of \f$ S\f$ is equal to the
orientation of the defining points, i.e., `orientation(A)`.

\cgalHeading{Implementation}

Spheres are implemented by a vector of points as a handle type. All
operations like creation, initialization, tests, input and output of a
sphere \f$ s\f$ take time \f$ O(s.dimension()) \f$. `dimension()`,
point access take constant time. The `center()`-operation takes
time \f$ O(d^3)\f$ on its first call and constant time thereafter. The
sidedness and orientation tests take time \f$ O(d^3)\f$. The space
requirement for spheres is \f$ O(s.dimension())\f$ neglecting the
storage room of the points.

*/
template< typename Kernel >
class Sphere_d {
public:

/// \name Types
/// @{

/*!
the linear algebra layer.
*/
typedef unspecified_type LA;

/*!
a read-only iterator for points defining
the sphere.
*/
typedef unspecified_type point_iterator;

/// @}

/// \name Creation
/// @{

/*!
introduces a variable `S`
of type `Sphere_d<Kernel>`.
*/
Sphere_d<Kernel>();

/*!
introduces a variable
`S` of type `Sphere_d<Kernel>`. `S` is initialized to the
sphere through the points in `A = tuple [first,last)`.

\pre \f$A\f$ consists of \f$d+1\f$ \f$d\f$-dimensional points.

\tparam ForwardIterator has `Point_d<Kernel>` as value type.
*/
template <class ForwardIterator> Sphere_d<Kernel>(int d,
ForwardIterator first, ForwardIterator last);

/// @}

/// \name Operations
/// @{

/*!
returns the dimension of the ambient space.
*/
int dimension();

/*!
returns the \f$ i\f$th defining
point.

\pre \f$ 0 \le i \le dim \f$
*/
Point_d<Kernel> point(int i) ;

/*!
returns an iterator
pointing to the first defining point.
*/
point_iterator points_begin() ;

/*!
returns an iterator pointing
beyond the last defining point.
*/
point_iterator points_end() ;

/*!
returns true iff the defining points
are not full dimensional.
*/
bool is_degenerate();

/*!
returns true iff the set of defining
points is legal. A set of defining points is legal iff their
orientation is non-zero or if they are all equal.
*/
bool is_legal() ;

/*!
returns the center of `S`.
\pre `S` is legal.
*/
Point_d<Kernel> center() ;

/*!
returns the squared radius of the
sphere.

\pre `S` is legal.
*/
FT squared_radius() ;

/*!
returns the orientation of
`S`.
*/
Orientation orientation() ;

/*!
returns
either the constant `ON_ORIENTED_BOUNDARY`,
`ON_POSITIVE_SIDE`, or `ON_NEGATIVE_SIDE`, iff p lies on the
boundary, properly on the positive side, or properly on the negative
side of sphere, resp.

\pre `S.dimension()==p.dimension()`.

*/
Oriented_side oriented_side(const Point_d<Kernel>& p) ;

/*!
returns
`ON_BOUNDED_SIDE`, `ON_BOUNDARY`, or `ON_UNBOUNDED_SIDE`
iff p lies properly inside, on the boundary, or properly outside of
sphere, resp.

\pre `S.dimension()==p.dimension()`.
*/
Bounded_side bounded_side(const Point_d<Kernel>& p) ;

/*!
returns
`S.oriented_side(p)==ON_POSITIVE_SIDE`.

\pre `S.dimension()==p.dimension()`.
*/
bool has_on_positive_side (const Point_d<Kernel>& p) ;

/*!
returns
`S.oriented_side(p)==ON_NEGATIVE_SIDE`.

\pre `S.dimension()==p.dimension()`.
*/
bool has_on_negative_side (const Point_d<Kernel>& p) ;

/*!
returns
`S.oriented_side(p)==ON_ORIENTED_BOUNDARY`, which is the same as
`S.bounded_side(p)==ON_BOUNDARY`.

\pre `S.dimension()==p.dimension()`.
*/
bool has_on_boundary (const Point_d<Kernel>& p) ;

/*!
returns
`S.bounded_side(p)==ON_BOUNDED_SIDE`.

\pre `S.dimension()==p.dimension()`.
*/
bool has_on_bounded_side (const Point_d<Kernel>& p) ;

/*!
returns
`S.bounded_side(p)==ON_UNBOUNDED_SIDE`.

\pre `S.dimension()==p.dimension()`.
*/
bool has_on_unbounded_side (const Point_d<Kernel>& p) ;

/*!
returns the sphere with the same
center and squared radius as `S` but with opposite orientation.

*/
Sphere_d<Kernel> opposite() ;

/*!
returns the
sphere translated by `v`.

\pre `S.dimension()==v.dimension()`.
*/
Sphere_d<Kernel> operator+(const Vector_d<Kernel>& v) ;

/// @}

}; /* end Sphere_d */

/*!
Test for equality as unoriented spheres.
\pre `S1.dimension()==S2.dimension()`.
\relates Sphere_d
*/
bool weak_equality(const Sphere_d<Kernel>& S1, const Sphere_d<Kernel>& S2) ;

} /* end namespace CGAL */

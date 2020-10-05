namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An instance of data type `Hyperplane_d` is an oriented hyperplane
in \f$ d\f$ - dimensional space. A hyperplane \f$ h\f$ is represented by
coefficients \f$ (c_0,c_1,\ldots,c_d)\f$ of type `RT`. At least one of
\f$ c_0\f$ to \f$ c_{ d - 1 }\f$ must be non-zero. The plane equation is
\f$ \sum_{ 0 \le i < d } c_i x_i + c_d = 0\f$, where \f$ x_0\f$ to \f$ x_{d-1}\f$ are
%Cartesian point coordinates. For a particular \f$ x\f$ the sign of \f$ \sum_{
0 \le i < d } c_i x_i + c_d\f$ determines the position of a point \f$ x\f$
with respect to the hyperplane (on the hyperplane, on the negative
side, or on the positive side).

There are two equality predicates for hyperplanes. The (weak) equality
predicate (`weak_equality`) declares two hyperplanes equal if they
consist of the same set of points, the strong equality predicate
(`operator==`) requires in addition that the negative halfspaces
agree. In other words, two hyperplanes are strongly equal if their
coefficient vectors are positive multiples of each other and they are
(weakly) equal if their coefficient vectors are multiples of each
other.

\cgalHeading{Implementation}

Hyperplanes are implemented by arrays of integers as an item type.
All operations like creation, initialization, tests, vector
arithmetic, input and output on a hyperplane \f$ h\f$ take time
\f$ O(h.dimension())\f$. coordinate access and `dimension()` take
constant time. The space requirement is \f$ O(h.dimension())\f$.

*/
template< typename Kernel >
class Hyperplane_d {
public:

/// \name Types
/// @{

/*!
the linear algebra layer.
*/
typedef unspecified_type LA;

/*!
a read-only iterator for the
coefficients.
*/
typedef unspecified_type Coefficient_const_iterator;

/// @}

/// \name Creation
/// @{

/*!
introduces a variable
`h` of type `Hyperplane_d<Kernel>`.
*/
Hyperplane_d<Kernel>();

/*!
introduces a
variable `h` of type `Hyperplane_d<Kernel>` initialized to the
hyperplane with coefficients `set [first,last)` and `D`.
\pre `size [first,last) == d`.
\tparam InputIterator has `RT` as value type.
*/
template <class InputIterator> Hyperplane_d<Kernel>(int d,
InputIterator first, InputIterator last, RT D);

/*!
introduces a variable
`h` of type `Hyperplane_d<Kernel>` initialized to the hyperplane
with coefficients `set [first,last)`.

\pre `size [first,last) == d+1`.
\tparam InputIterator has `RT` as value type.
*/
template <class InputIterator> Hyperplane_d<Kernel>(int d,
InputIterator first, InputIterator last);

/*!
constructs
some hyperplane that passes through the points in `set [first,last)`. If `side` is `ON_POSITIVE_SIDE` or
`ON_NEGATIVE_SIDE` then `o` is on that side of the
constructed hyperplane.

\pre A hyperplane with the stated properties must exist.
\tparam ForwardIterator has `Point_d<Kernel>` as value type.
*/
template <class ForwardIterator>
Hyperplane_d<Kernel>(ForwardIterator first, ForwardIterator last,
Point_d<Kernel> o, Oriented_side side = ON_ORIENTED_BOUNDARY);

/*!
constructs the hyperplane with normal direction `dir` that
passes through \f$ p\f$. The direction `dir` points into the positive
side.

\pre `p.dimension()==dir.dimension()` and `dir` is not degenerate.
*/
Hyperplane_d<Kernel>(Point_d<Kernel> p, Direction_d<Kernel>
dir);

/*!
introduces a
variable `h` of type `Hyperplane_d<Kernel>` in \f$ 2\f$-dimensional
space with equation \f$ ax+by+c=0\f$.
*/
Hyperplane_d<Kernel>(RT a, RT b, RT c);

/*!
introduces a
variable `h` of type `Hyperplane_d<Kernel>` in \f$ 3\f$-dimensional
space with equation \f$ ax+by+cz+d=0\f$.
*/
Hyperplane_d<Kernel>(RT a, RT b, RT c, RT d);

/// @}

/// \name Operations
/// @{

/*!
returns the dimension of `h`.
*/
int dimension() ;

/*!
returns the \f$ i\f$-th coefficient of
`h`.

\pre \f$ 0 \leq i \leq d\f$.
*/
RT operator[](int i) ;

/*!
returns the \f$ i\f$-th coefficient of
`h`.

\pre \f$ 0 \leq i \leq d\f$.
*/
RT coefficient(int i) ;

/*!
returns
an iterator pointing to the first coefficient.
*/
Coefficient_const_iterator coefficients_begin() ;

/*!
returns an
iterator pointing beyond the last coefficient.
*/
Coefficient_const_iterator coefficients_end() ;

/*!
returns the orthogonal
vector of `h`. It points from the negative halfspace into the
positive halfspace and its homogeneous coordinates are \f$ (c_0,
\ldots, c_{d - 1},1)\f$.
*/
Vector_d<Kernel> orthogonal_vector() ;

/*!
returns the
orthogonal direction of `h`. It points from the negative
halfspace into the positive halfspace.
*/
Direction_d<Kernel> orthogonal_direction() ;

/*!
returns
the side of the hyperplane `h` containing \f$ p\f$.

\pre `h.dimension() == p.dimension()`.
*/
Oriented_side oriented_side(const Point_d<Kernel>& p) ;

/*!
returns true iff point
`p` lies on the hyperplane `h`.

\pre `h.dimension() == p.dimension()`.
*/
bool has_on(const Point_d<Kernel>& p) ;

/*!
returns true
iff point `p` lies on the boundary of hyperplane `h`.
\pre `h.dimension() == p.dimension()`.
*/
bool has_on_boundary(const Point_d<Kernel>& p) ;

/*!
returns
true iff point `p` lies on the positive side of hyperplane
`h`.

\pre `h.dimension() == p.dimension()`.
*/
bool has_on_positive_side(const Point_d<Kernel>& p) ;

/*!
returns
true iff point `p` lies on the negative side of hyperplane
`h`.

\pre `h.dimension() == p.dimension()`.
*/
bool has_on_negative_side(const Point_d<Kernel>& p) ;

/*!
returns \f$ t(h)\f$.

\pre `h.dimension() == t.dimension()`.

*/
Hyperplane_d<Kernel> transform(const Aff_transformation_d<Kernel>& t)
;

/// @}

}; /* end Hyperplane_d */

/*!
test for weak equality.

\pre `h1.dimension() == h2.dimension()`.
\relates Hyperplane_d
*/
bool weak_equality(const Hyperplane_d<Kernel>& h1, const Hyperplane_d<Kernel>& h2) ;

} /* end namespace CGAL */

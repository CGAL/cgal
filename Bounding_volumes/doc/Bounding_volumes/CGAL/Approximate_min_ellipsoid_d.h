
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

An object of class `Approximate_min_ellipsoid_d` is an approximation to the
ellipsoid of smallest volume enclosing a finite multiset of points
in \f$ d\f$-dimensional Euclidean space \f$ \E^d\f$, \f$ d\ge 2\f$.

An <I>ellipsoid</I> in \f$ \E^d\f$ is a Cartesian pointset of the form \f$ \{
x\in\E^d \mid x^T E x + x^T e + \eta\leq 0 \}\f$, where \f$ E\f$ is some
positive definite matrix from the set \f$ \mathbb{R}^{d\times d}\f$, \f$ e\f$ is some
real \f$ d\f$-vector, and \f$ \eta\in\mathbb{R}\f$. A pointset \f$ P\subseteq \E^d\f$ is
called <I>full-dimensional</I> if its affine hull has dimension \f$ d\f$.
For a finite, full-dimensional pointset \f$ P\f$ we denote by \f$ (P)\f$ the
smallest ellipsoid that contains all points of \f$ P\f$; this ellipsoid
exists and is unique.

For a given finite and full-dimensional pointset \f$ P\subset \E^d\f$ and a
real number \f$ \epsilon\ge 0\f$, we say that an ellipsoid \f$ {\cal
E}\subset\E^d\f$ is an <I>\f$ (1+\epsilon)\f$-appoximation</I> to \f$ (P)\f$ if
\f$ P\subset {\cal E}\f$ and \f$ ({\cal E}) \leq (1+\epsilon)
((P))\f$. In other words, an \f$ (1+\epsilon)\f$-approximation to
\f$ (P)\f$ is an enclosing ellipsoid whose volume is by at most a
factor of \f$ 1+\epsilon\f$ larger than the volume of the smallest
enclosing ellipsoid of \f$ P\f$.

Given this notation, an object of class `Approximate_min_ellipsoid_d` represents an
\f$ (1+\epsilon)\f$-approximation to \f$ (P)\f$ for a given finite and
full-dimensional multiset of points \f$ P\subset\E^d\f$ and a real constant
\f$ \epsilon>0\f$.\cgalFootnote{A <I>multiset</I> is a set where elements may have multiplicity greater than \f$ 1\f$.} When an
`Approximate_min_ellipsoid_d<Traits>` object is constructed, an
iterator over the points \f$ P\f$ and the number \f$ \epsilon\f$ have to be
specified; the number \f$ \epsilon\f$ defines the <I>desired
approximation ratio</I> \f$ 1+\epsilon\f$. The underlying algorithm will then
try to compute an \f$ (1+\epsilon)\f$-approximation to \f$ (P)\f$, and one of
the following two cases takes place.
<UL>
<LI>The algorithm determines that \f$ P\f$ is not full-dimensional (see
`is_full_dimensional()` below).

<I>Important note:</I> due to rounding errors, the algorithm cannot
in all cases decide correctly whether \f$ P\f$ is full-dimensional or
not. If `is_full_dimensional()` returns `false`, the points
lie in such a "thin" subspace of \f$ \E^d\f$ that the algorithm is
incapable of computing an approximation to \f$ (P)\f$. More
precisely, if `is_full_dimensional()` returns `false`, there
exist two parallel hyperplanes in \f$ \E^d\f$ with the points \f$ P\f$ in
between so that the distance \f$ \delta\f$ between the hyperplanes is
very small, possible zero. (If \f$ \delta=0\f$ then \f$ P\f$ is not
full-dimensional.)

If \f$ P\f$ is not full-dimensional, linear algebra techniques should be
used to determine an affine subspace \f$ S\f$ of \f$ \E^d\f$ that contains the
points \f$ P\f$ as a (w.r.t.\ \f$ S\f$) full-dimensional pointset; once \f$ S\f$ is
determined, the algorithm can be invoked again to compute an
approximation to (the lower-dimensional) \f$ (P)\f$ in \f$ S\f$. Since
`is_full_dimensional()` might (due to rounding errors, see
above) return `false` even though \f$ P\f$ is full-dimensional, the
lower-dimensional subspace \f$ S\f$ containing \f$ P\f$ need not exist.
Therefore, it might be more advisable to fit a hyperplane \f$ H\f$
through the pointset \f$ P\f$, project \f$ P\f$ onto this affine subspace \f$ H\f$,
and compute an approximation to the minimum-volume enclosing
ellipsoid of the projected points within \f$ H\f$; the fitting can be
done for instance using the `linear_least_squares_fitting()`
function from the \cgal package `Principal_component_analysis`.
<LI>The algorithm determines that \f$ P\f$ is full-dimensional. In this
case, it provides an approximation \f$ {\cal E}\f$ to \f$ (P)\f$, but
depending on the input problem (i.e., on the pair \f$ (P,\epsilon)\f$),
it may not have achieved the desired approximation ratio but merely
some <I>worse</I> approximation ratio \f$ 1+\epsilon'>1+\epsilon\f$. The
achieved approximation ratio \f$ 1+\epsilon'\f$ can be queried using
`achieved_epsilon()`, which returns \f$ \epsilon'\f$. The ellipsoid
\f$ {\cal E}\f$ itself can be queried via the methods
`defining_matrix()`, `defining_vector()`, and
`defining_scalar()`.
</UL>

The ellipsoid \f$ {\cal E}\f$ computed by the algorithm satisfies the inclusions
\anchor eqapproximate_min_ellipsoid_incl
\f[
\frac{1}{(1+\epsilon')d} {\cal E} \subseteq \mathop{\rm conv}\nolimits(P) \subseteq {\cal E}
\f]
where \f$ f {\cal E}\f$ denotes the ellipsoid \f$ {\cal E}\f$ scaled by the
factor \f$ f\in\mathbb{R}^+\f$ with respect to its center, and where \f$ \mathop{\rm
conv}\nolimits(A)\f$ denotes the <I>convex hull</I> of a pointset
\f$ A\subset \E^d\f$.

The underlying algorithm can cope with all kinds of inputs (multisets
\f$ P\f$, \f$ \epsilon\in[0,\infty)\f$) and terminates in all cases. There
is, however, no guarantee that any desired approximation ratio
is actually achieved; the performance of the algorithm in this respect
highly depends on the input pointset. Values of at least \f$ 0.01\f$ for
\f$ \epsilon\f$ are usually handled without problems.

Internally, the algorithm represents the input points' Cartesian
coordinates as `double`'s. For this conversion to work, the input
point coordinates must be convertible to `double`. Also, in order
to compute the achieved epsilon \f$ \epsilon'\f$ mentioned above, the algorithm
requires a number type `ET` that provides <I>exact</I> arithmetic.
(Both these aspects are discussed in the documentation of the concept
`ApproximateMinEllipsoid_d_Traits_d`.)


\tparam Traits must be a model for `ApproximateMinEllipsoid_d_Traits_d`.

We provide the model `CGAL::Approximate_min_ellipsoid_d_traits_d<K>`
using the \f$ d\f$-dimensional \cgal kernel; the models
`CGAL::Approximate_min_ellipsoid_d_traits_2<K>` and
`CGAL::Approximate_min_ellipsoid_d_traits_3<K>` are for use with the
\f$ 2\f$- and \f$ 3\f$-dimensional \cgal kernel, respectively.

\sa `CGAL::Min_ellipse_2<Traits>`

\cgalHeading{Implementation}

We implement Khachyian's algorithm for rounding
polytopes \cgalCite{cgal:k-rprnm-96}. Internally, we use
`double`-arithmetic and (initially a single)
Cholesky-decomposition. The algorithm's running time is
\f$ {\cal O}(nd^2(\epsilon^{-1}+\ln d + \ln\ln(n)))\f$, where \f$ n=|P|\f$ and
\f$ 1+\epsilon\f$ is the desired approximation ratio.

\cgalHeading{Example}

To illustrate the usage of `Approximate_min_ellipsoid_d` we give two examples in 2D. The
first program generates a random set \f$ P\subset\E^2\f$ and outputs the
points and a \f$ 1.01\f$-approximation of \f$ (P)\f$ as an EPS-file, which
you can view using <TT>gv</TT>, for instance. (In both examples you can
change the variables `n` and `d` to experiment with the code.)

\cgalExample{Approximate_min_ellipsoid_d/ellipsoid.cpp}

The second program outputs the approximation in a format suitable
for display in Maplesoft's Maple.

\cgalExample{Approximate_min_ellipsoid_d/ellipsoid_for_maple.cpp}

\note This class requires the \ref thirdpartyEigen library.

*/
template< typename Traits >
class Approximate_min_ellipsoid_d {
public:

/// \name Types
/// @{

/*!
`typedef Traits::FT FT` (which is always a
typedef to `double`).
*/
typedef unspecified_type FT;

/*!
`typedef Traits::ET ET` (which is an exact number type used for exact computation like for example in `achieved_epsilon()`).
*/
typedef unspecified_type ET;

/*!
`typedef Traits::Point Point`
*/
typedef unspecified_type Point;

/*!
`typedef Traits::Cartesian_const_iterator Cartesian_const_iterator`
*/
typedef unspecified_type Cartesian_const_iterator;

/*!
A model of STL concept
`RandomAccessIterator` with value type `double` that is used
to iterate over the Cartesian center coordinates of the computed
ellipsoid, see `center_cartesian_begin()`.
*/
typedef unspecified_type Center_coordinate_iterator;

/*!
A model of STL concept
`RandomAccessIterator` with value type `double` that is used
to iterate over the lengths of the semiaxes of the computed ellipsoid,
see `axes_lengths_begin()`.
*/
typedef unspecified_type Axes_lengths_iterator;

/*!
A model of STL concept
`RandomAccessIterator` with value type `double` that is used
to iterate over the %Cartesian coordinates of the direction of a fixed
axis of the computed ellipsoid, see
`axis_direction_cartesian_begin()`.
*/
typedef unspecified_type Axes_direction_coordinate_iterator;

/// @}

/// \name Creation
/// An object of type `Approximate_min_ellipsoid_d` can be created
/// from an arbitrary point set \f$ P\f$ and some nonnegative `double`
/// value `eps`.
/// @{

/*!

initializes `ame` to an \f$ (1+\epsilon)\f$-approximation of
\f$ (P)\f$ with \f$ P\f$ being the set of points in the range
[`first`,`last`). The number \f$ \epsilon\f$ in this will
be at most `eps`, if possible. However, due to the
limited precision in the algorithm's underlying arithmetic, it
can happen that the computed approximation ellipsoid has a
worse approximation ratio (and \f$ \epsilon\f$ can thus be larger
than `eps` in general). In any case, the number
\f$ \epsilon\f$ (and with this, the achived approximation
\f$ 1+\epsilon\f$) can be queried by calling the routine
`achieved_epsilon()` discussed below.

\tparam Iterator must be a model of `InputIterator` with `Point` as value type.

\pre The dimension \f$ d\f$ of the input points must be at least \f$ 2\f$, and \f$ \epsilon>0\f$.
*/
template < class Iterator >
Approximate_min_ellipsoid_d(double eps,
Iterator first,
Iterator last,
const Traits& traits = Traits() );

/// @}

/// \name Access Functions
/// The following methods can be used to query the achieved
/// approximation ratio \f$ 1+\epsilon'\f$ and the computed ellipsoid
/// \f$ {\cal E} = \{ x\in\E^d \mid x^T E x + x^T e + \eta\leq 0
/// \}\f$. The methods `defining_matrix()`, `defining_vector()`, and
/// `defining_scalar()` do not return \f$ E\f$, \f$ e\f$, and \f$
/// \eta\f$ directly but yield multiples of these quantities that are
/// exactly representable using the `double` type. (This is necessary
/// because the parameters \f$ E\f$, \f$ e\f$, and \f$ \eta\f$ of the
/// computed approximation ellipsoid \f$ {\cal E}\f$ might not be
/// exactly representable as `double` numbers.) In order to access the
/// center and semiaxes of the computed approximation ellipsoid, the
/// functions `center_cartesian_begin()`, `axes_lengths_begin()`, and
/// `axis_direction_cartesian_begin()` can be used. In constrast to
/// the above access functions `achieved_epsilon()`,
/// `defining_matrix()`, `defining_vector()`, and `defining_scalar()`,
/// which return the described quantities exactly, the routines below
/// return <I>numerical approximations</I> to the real center and real
/// semiaxes of the computed ellipsoid; the comprised relative error
/// may be larger than zero, and there are no guarantees for the
/// returned quantities.
/// @{

/*!

returns the number of points of `ame`, i.e., \f$ |P|\f$.
*/
unsigned int number_of_points( ) const;

/*!
returns a number
\f$ \epsilon'\f$ such that the computed approximation is (under exact
arithmetic) guaranteed to be an \f$ (1+\epsilon')\f$-approximation to
\f$ (P)\f$.
\pre `ame.is_full_dimensional() == true`.
\post \f$ \epsilon'>0\f$.
*/
double achieved_epsilon() const;

/*!

gives access to the \f$ (i,j)\f$th entry of the matrix \f$ E\f$ in the
representation \f$ \{ x\in\E^d \mid x^T E x + x^T e + \eta\leq0
\}\f$ of the computed approximation ellipsoid \f$ {\cal E}\f$. The number returned by
this routine is \f$ (1+\epsilon')(d+1)\,E_{ij}\f$, where \f$ \epsilon'\f$ is
the number returned by `achieved_epsilon()`.
\pre \f$ 0\leq i,j\leq d\f$, where \f$ d\f$ is the dimension of the points \f$ P\f$, and `ame.is_full_dimensional() == true`.
*/
double defining_matrix(int i,int j) const;

/*!

gives access to the \f$ i\f$th entry of the vector \f$ e\f$ in the
representation \f$ \{ x\in\E^d \mid x^T E x + x^T e + \eta\leq0
\}\f$ of the computed approximation ellipsoid \f$ {\cal E}\f$. The number returned by
this routine is \f$ (1+\epsilon')(d+1)\,e_{i}\f$, where \f$ \epsilon'\f$ is
the number returned by `achieved_epsilon()`.
\pre \f$ 0\leq i\leq d\f$, where \f$ d\f$ is the dimension of the points \f$ P\f$, and `ame.is_full_dimensional() == true`.
*/
double defining_vector(int i) const;

/*!

gives access to the scalar \f$ \eta\f$ from the
representation \f$ \{ x\in\E^d \mid x^T E x + x^T e + \eta\leq0
\}\f$ of the computed approximation ellipsoid \f$ {\cal E}\f$. The number returned by
this routine is \f$ (1+\epsilon')(d+1)\,(\eta+1)\f$, where \f$ \epsilon'\f$ is
the number returned by `achieved_epsilon()`.
\pre `ame.is_full_dimensional() == true`.
*/
double defining_scalar() const;

/*!

returns a const reference to the traits class object.
*/
const Traits& traits( ) const;

/*!
returns the dimension of
the ambient space, i.e., the dimension of the points \f$ P\f$.
*/
int dimension() const;

/*!

returns an iterator pointing to the first of the \f$ d\f$ Cartesian
coordinates of the computed ellipsoid's center.

The returned point is a floating-point approximation to the
ellipsoid's exact center; no guarantee is given w.r.t.\ the involved
relative error. \pre `ame.is_full_dimensional() == true`.
*/
Center_coordinate_iterator center_cartesian_begin();

/*!

returns the past-the-end iterator corresponding to
`center_cartesian_begin()`.
\pre `ame.is_full_dimensional() == true`.
*/
Center_coordinate_iterator center_cartesian_end();

/*!

returns an iterator pointing to the first of the \f$ d\f$ descendantly
sorted lengths of the computed ellipsoid's axes. The \f$ d\f$ returned
numbers are floating-point approximations to the exact
axes-lengths of the computed ellipsoid; no guarantee is given
w.r.t.\ the involved relative error. (See also method
`axes_direction_cartesian_begin()`.) \pre `ame.is_full_dimensional() == true`, and \f$ d\in\{2,3\}\f$.
*/
Axes_lengths_iterator axes_lengths_begin();

/*!

returns the past-the-end iterator corresponding to
`axes_lengths_begin()`. \pre `ame.is_full_dimensional() == true`, and \f$ d\in\{2,3\}\f$.
*/
Axes_lengths_iterator axes_lengths_end();

/*!
returns an iterator
pointing to the first of the \f$ d\f$ %Cartesian coordinates of the
computed ellipsoid's \f$ i\f$th axis direction (i.e., unit vector in
direction of the ellipsoid's \f$ i\f$th axis). The direction described
by this iterator is a floating-point approximation to the exact
axis direction of the computed ellipsoid; no guarantee is given
w.r.t.\ the involved relative error. An approximation to the
length of axis \f$ i\f$ is given by the \f$ i\f$th entry of
`axes_lengths_begin()`.
\pre `ame.is_full_dimensional() == true`, and \f$ d\in\{2,3\}\f$, and \f$ 0\leq i < d\f$.
*/
Axes_direction_coordinate_iterator
axis_direction_cartesian_begin(int i);

/*!
returns the past-the-end
iterator corresponding to
`axis_direction_cartesian_begin()`.\pre `ame.is_full_dimensional() == true`, and \f$ d\in\{2,3\}\f$, and \f$ 0\leq i < d\f$.
*/
Axes_direction_coordinate_iterator
axis_direction_cartesian_end(int i);

/// @}

/// \name Predicates
/// @{

/*!
returns whether
\f$ P\f$ is full-dimensional or not, i.e., returns `true` if
and only if \f$ P\f$ is full-dimensional.
<I>Note:</I> due to
the limited precision in the algorithm's underlying
arithmetic, the result of this method is not always
correct. Rather, a return value of `false` means that the
points \f$ P\f$ are contained in a "very thin" linear subspace of
\f$ \E^d\f$, and as a consequence, the algorithm cannot compute an
approximation. More precisely, a return value of `false`
means that the points \f$ P\f$ are contained between two parallel
hyperplanes in \f$ \E^d\f$ that are very close to each other
(possibly at distance zero) - so close, that the algorithm
could not compute an approximation ellipsoid. Similarly, a
return value of `true` does not guarantee \f$ P\f$ to be
full-dimensional; but there exists an input pointset \f$ P'\f$ such
that the points \f$ P'\f$ and \f$ P\f$ have almost identical coordinates
and \f$ P'\f$ is full-dimensional.
*/
bool is_full_dimensional( ) const;

/// @}

/// \name Validity Check
/// An object `ame` is valid iff <UL> <LI>`ame` contains all points of
/// its defining set \f$ P\f$, <LI>`ame` is an \f$
/// (1+\epsilon')\f$-approximation to the smallest ellipsoid \f$
/// (P)\f$ of \f$ P\f$, <LI>The ellipsoid represented by `ame`
/// fulfills the inclusion ( \ref eqapproximate_min_ellipsoid_incl
/// ). </UL>
/// @{

/*!

returns `true` iff `ame` is valid according to the above
definition. If `verbose` is `true`, some messages
concerning the performed checks are written to the standard error
stream.
*/
bool is_valid( bool verbose = false) const;

/// @}

/// \name Miscellaneous
/// @{

/*!

Writes the points \f$ P\f$ and the computed approximation to
\f$ (P)\f$ as an EPS-file under pathname `name`. \pre The dimension of points \f$ P\f$ must be \f$ 2\f$.
<I>Note:</I> this
routine is provided as a debugging routine; future version of
\cgal might not provide it anymore.
\pre `ame.is_full_dimensional() == true`.
*/
void write_eps(const std::string& name) const;

/// @}

}; /* end Approximate_min_ellipsoid_d */
} /* end namespace CGAL */

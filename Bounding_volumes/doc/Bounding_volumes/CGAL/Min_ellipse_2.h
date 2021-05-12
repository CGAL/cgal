
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

An object of the class `Min_ellipse_2` is the unique ellipse of smallest area
enclosing a finite (multi)set of points in two-dimensional euclidean
space \f$ \E^2\f$. For a point set \f$ P\f$ we denote by \f$ me(P)\f$ the smallest
ellipse that contains all points of \f$ P\f$. Note that \f$ me(P)\f$ can be
degenerate,
i.e.\ \f$ me(P)=\emptyset\f$ if
\f$ P=\emptyset\f$, \f$ me(P)=\{p\}\f$ if \f$ P=\{p\}\f$,
and <span class="mbox">\f$ me(P) = \{ (1-\lambda)p + \lambda q \mid 0 \leq \lambda \leq 1 \}\f$</span> if \f$ P=\{p,q\}\f$.

An inclusion-minimal subset \f$ S\f$ of \f$ P\f$ with \f$ me(S)=me(P)\f$ is called a
<I>support set</I>,
the points in \f$ S\f$ are the <I>support points</I>. A support set has size at
most five, and all its points lie on the boundary of \f$ me(P)\f$. In general,
neither the support set nor its size are necessarily unique.

The underlying algorithm can cope with all kinds of input, e.g. \f$ P\f$ may be
empty or points may occur more than once. The algorithm computes a support
set \f$ S\f$ which remains fixed until the next insert or clear operation.

\tparam Traits must be a model for `MinEllipse2Traits`.

We provide the model `CGAL::Min_ellipse_2_traits_2<K>` using the
two-dimensional \cgal kernel.

\sa `CGAL::Min_circle_2<Traits>`
\sa `CGAL::Min_ellipse_2_traits_2<K>`
\sa `MinEllipse2Traits`

\cgalHeading{Implementation}

We implement the incremental algorithm of Welzl, with move-to-front
heuristic \cgalCite{w-sedbe-91a}, using the primitives as described
in \cgalCite{gs-epsee-97}, \cgalCite{cgal:gs-seefe-97a}. The whole implementation is described
in \cgalCite{cgal:gs-seeeg-98}.

If randomization is
chosen, the creation time is almost always linear in the number of points.
Access functions and predicates take constant time, inserting a point might
take up to linear time, but substantially less than computing the new
smallest enclosing ellipse from scratch. The clear operation and the check
for validity each takes linear time.

\cgalHeading{Example}

To illustrate the usage of `Min_ellipse_2` and to show that randomization
can be useful in certain cases, we give an example. The example also
shows how the coefficents of the constructed ellipse can be accessed.

\cgalExample{Min_ellipse_2/min_ellipse_2.cpp}

*/
template< typename Traits >
class Min_ellipse_2 {
public:

/// \name Types
/// @{

/*!
Typedef to `Traits::Point `.
*/
typedef unspecified_type Point ;

/*!
Typedef to `Traits::Ellipse`. If you
are using the predefined traits class
`CGAL::Min_ellipse_2_traits_2<K>`,
you can access the coefficients of the ellipse, see the
documentation of `CGAL::Min_ellipse_2_traits_2<K>` and
the example below.
*/
typedef unspecified_type Ellipse;

/*!

Non-mutable model of the \stl concept <I>BidirectionalIterator</I>
with value type `Point`. Used to access the points
of the smallest enclosing ellipse.
*/
typedef unspecified_type Point_iterator;

/*!

Non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `Point`. Used to access the support points
of the smallest enclosing ellipse.
*/
typedef unspecified_type Support_point_iterator;

/// @}

/// \name Creation
/// A `Min_ellipse_2` object can be created from an arbitrary point
/// set \f$ P\f$ and by specialized construction methods expecting no,
/// one, two, three, four or five points as arguments. The latter
/// methods can be useful for reconstructing \f$ me(P)\f$ from a given
/// support set \f$ S\f$ of \f$ P\f$.
/// @{

/*!

initializes `min_ellipse` to \f$ me(P)\f$ with \f$ P\f$ being the set of points
in the range [`first`,`last`). If `randomize` is
`true`, a random permutation of \f$ P\f$ is computed in
advance, using the random numbers generator `random`.
Usually, this will not be necessary, however, the algorithm's
efficiency depends on the order in which the points are
processed, and a bad order might lead to extremely poor
performance (see example below).
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
*/
template < class InputIterator >
Min_Ellipse_2( InputIterator first,
InputIterator last,
bool randomize,
Random& random = get_default_random(),
const Traits& traits = Traits() );

/*!

creates a variable `min_ellipse` of type `Min_ellipse_2`.
It is initialized to
\f$ me(\emptyset)\f$, the empty set.
\post `min_ellipse.is_empty()` = `true`.
*/
Min_ellipse_2( const Traits& traits = Traits());

/*!

initializes `min_ellipse` to \f$ me(\{p\})\f$, the set \f$ \{p\}\f$.
\post `min_ellipse.is_degenerate()` = `true`.
*/
Min_ellipse_2( const Point& p,
const Traits& traits = Traits());

/*!

initializes `min_ellipse` to \f$ me(\{p,q\})\f$,
the set \f$ \{ (1-\lambda) p + \lambda q \mid 0 \leq \lambda \leq 1 \}\f$
\post `min_ellipse.is_degenerate()` = `true`.
*/
Min_ellipse_2( const Point& p,
const Point& q,
const Traits& traits = Traits());

/*!

initializes `min_ellipse` to \f$ me(\{p1,p2,p3\})\f$.
*/
Min_ellipse_2( const Point& p1,
const Point& p2,
const Point& p3,
const Traits& traits = Traits());

/*!

initializes `min_ellipse` to \f$ me(\{p1,p2,p3,p4\})\f$.
*/
Min_ellipse_2( const Point& p1,
const Point& p2,
const Point& p3,
const Point& p4,
const Traits& traits = Traits());

/*!

initializes `min_ellipse` to \f$ me(\{p1,p2,p3,p4,p5\})\f$.
*/
Min_ellipse_2( const Point& p1,
const Point& p2,
const Point& p3,
const Point& p4,
const Point& p5,
const Traits& traits = Traits());

/// @}

/// \name Access Functions
/// @{

/*!

returns the number of points of `min_ellipse`, i.e.\ \f$ |P|\f$.
*/
int number_of_points( ) const;

/*!

returns the number of support points of `min_ellipse`, i.e.\ \f$ |S|\f$.
*/
int number_of_support_points( ) const;

/*!

returns an iterator referring to the first point of `min_ellipse`.
*/
Point_iterator points_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Point_iterator points_end() const;

/*!

returns an iterator referring to the first support point of `min_ellipse`.
*/
Support_point_iterator support_points_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Support_point_iterator support_points_end() const;

/*!

returns the `i`-th support point of `min_ellipse`. Between two
modifying operations (see below) any call to
`min_ellipse.support_point(i)` with the same `i` returns
the same point.
\pre \f$ 0 \leq i<\f$ `min_ellipse.number_of_support_points()`.
*/
const Point& support_point( int i) const;

/*!

returns the current ellipse of `min_ellipse`.
*/
const Ellipse& ellipse( ) const;

/// @}

/// \name Predicates
/// By definition, an empty `Min_ellipse_2` has no boundary and no
/// bounded side, i.e.\ its unbounded side equals the whole space \f$
/// \E^2\f$.
/// @{

/*!

returns `CGAL::ON_BOUNDED_SIDE`,
`CGAL::ON_BOUNDARY`, or
`CGAL::ON_UNBOUNDED_SIDE` iff `p` lies properly
inside, on the boundary of, or properly outside of `min_ellipse`, resp.
*/
CGAL::Bounded_side
bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies properly inside `min_ellipse`.
*/
bool has_on_bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies on the boundary
of `min_ellipse`.
*/
bool has_on_boundary( const Point& p) const;

/*!

returns `true`, iff `p` lies properly outside of `min_ellipse`.
*/
bool has_on_unbounded_side( const Point& p) const;

/*!

returns `true`, iff `min_ellipse` is empty (this implies
degeneracy).
*/
bool is_empty( ) const;

/*!

returns `true`, iff `min_ellipse` is degenerate,
i.e.\ if `min_ellipse` is empty, equal to a single point or equal to a
segment, equivalently if the number of support points is less
than 3.
*/
bool is_degenerate( ) const;

/// @}

/// \name Modifiers
/// New points can be added to an existing `min_ellipse`, allowing to
/// build \f$ me(P)\f$ incrementally, e.g. if \f$ P\f$ is not known in
/// advance. Compared to the direct creation of \f$ me(P)\f$, this is
/// not much slower, because the construction method is incremental
/// itself.
/// @{

/*!

inserts `p` into `min_ellipse` and recomputes the smallest
enclosing ellipse.
*/
void insert( const Point& p);

/*!

inserts the points in the range [`first`,`last`)
into `min_ellipse` and recomputes the smallest enclosing ellipse by
calling `insert(p)` for each point `p` in
[`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
*/
template < class InputIterator >
void insert( InputIterator first,
InputIterator last );

/*!

deletes all points in `min_ellipse` and sets `min_ellipse` to the empty set.
\post `min_ellipse.is_empty()` = `true`.
*/
void clear( );

/// @}

/// \name Validity Check
/// An object `min_ellipse` is valid, iff <UL> <LI>`min_ellipse`
/// contains all points of its defining set \f$ P\f$,
/// <LI>`min_ellipse` is the smallest ellipse spanned by its support
/// set \f$ S\f$, and <LI>\f$ S\f$ is minimal, i.e.\ no support point
/// is redundant. </UL> <I>Note:</I> In this release only the first
/// item is considered by the validity check.
/// @{

/*!

returns `true`, iff `min_ellipse` contains all points of its
defining set \f$ P\f$. If `verbose` is `true`, some
messages concerning the performed checks are written to
standard error stream. The second parameter `level` is not
used, we provide it only for consistency with interfaces of
other classes.
*/
bool is_valid( bool verbose = false,
int level = 0 ) const;

/// @}

/// \name Miscellaneous
/// @{

/*!

returns a const reference to the traits class object.
*/
const Traits& traits( ) const;


/// @}

}; /* end Min_ellipse_2 */


/*!

writes `min_ellipse` to output stream `os`.
An overload of `operator<<` must be defined for `Point` (and for `Ellipse`, if pretty printing is used).
\relates Min_ellipse_2
*/
std::ostream&
operator << ( std::ostream& os,
const Min_ellipse_2<Traits>& min_ellipse);

/*!

reads `min_ellipse` from input stream `is`.
An overload of `operator>>` must be defined for `Point`.
\relates Min_ellipse_2
*/
std::istream&
operator >> ( std::istream& is,
Min_ellipse_2<Traits> min_ellipse&);

} /* end namespace CGAL */


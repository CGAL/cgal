
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

An object of the class `Min_circle_2` is the unique circle of smallest area
enclosing a finite (multi)set of points in two-dimensional Euclidean
space \f$ \E^2\f$. For a point set \f$ P\f$ we denote by \f$ mc(P)\f$ the smallest circle
that contains all points of \f$ P\f$. Note that \f$ mc(P)\f$ can be
degenerate,
i.e.\ \f$ mc(P)=\emptyset\f$ if
\f$ P=\emptyset\f$ and \f$ mc(P)=\{p\}\f$ if
\f$ P=\{p\}\f$.

An inclusion-minimal subset \f$ S\f$ of \f$ P\f$ with \f$ mc(S)=mc(P)\f$ is called a
<I>support set</I>,
the points in \f$ S\f$ are the <I>support points</I>. A support set has size at
most three, and all its points lie on the boundary of \f$ mc(P)\f$. In general,
neither the support set nor its size are necessarily unique.

The underlying algorithm can cope with all kinds of input, e.g. \f$ P\f$ may be
empty or points may occur more than once. The algorithm computes a support
set \f$ S\f$ which remains fixed until the next insert or clear operation.

\note This class is (almost) obsolete. The class
`CGAL::Min_sphere_of_spheres_d<Traits>` solves a more general problem
and is faster than `Min_circle_2` even if used only for points in two
dimensions as input. Most importantly,
`CGAL::Min_sphere_of_spheres_d<Traits>` has
a specialized implementation for floating-point arithmetic which
ensures correct results in a large number of cases (including
highly degenerate ones). In contrast, `Min_circle_2` is not tuned for
floating-point computations. The only advantage of
`Min_circle_2` over `CGAL::Min_sphere_of_spheres_d<Traits>` is that the
former can deal with points in homogeneous coordinates, in which
case the algorithm is division-free. Thus, `Min_circle_2` might still
be an option in case your input number type cannot (efficiently)
divide.


\tparam Traits must be a model for `MinCircle2Traits`.

We provide the model `CGAL::Min_circle_2_traits_2` using the
two-dimensional \cgal kernel.

\sa `CGAL::Min_ellipse_2<Traits>`
\sa `CGAL::Min_sphere_d<Traits>`
\sa `CGAL::Min_sphere_of_spheres_d<Traits>`
\sa `CGAL::Min_circle_2_traits_2<K>`
\sa `MinCircle2Traits`

\cgalHeading{Implementation}

We implement the incremental algorithm of Welzl, with move-to-front
heuristic \cgalCite{w-sedbe-91a}. The whole implementation is described
in \cgalCite{cgal:gs-seceg-98}.

If randomization is
chosen, the creation time is almost always linear in the number of points.
Access functions and predicates take constant time, inserting a point might
take up to linear time, but substantially less than computing the new
smallest enclosing circle from scratch. The clear operation and the check
for validity each takes linear time.

\cgalHeading{Example}

To illustrate the creation of `Min_circle_2` and to show that
randomization can be useful in certain cases, we give an example.

\cgalExample{Min_circle_2/min_circle_2.cpp}

*/
template< typename Traits >
class Min_circle_2 {
public:

/// \name Types
/// @{

/*!
typedef to `Traits::Point `.
*/
typedef unspecified_type Point ;

/*!
typedef to `Traits::Circle`.
*/
typedef unspecified_type Circle;

/*!

non-mutable model of the \stl concept <I>BidirectionalIterator</I>
with value type `Point`. Used to access the points
of the smallest enclosing circle.
*/
typedef unspecified_type Point_iterator;

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `Point`. Used to access the support points
of the smallest enclosing circle.
*/
typedef unspecified_type Support_point_iterator;

/// @}

/// \name Creation
/// A `Min_circle_2` object can be created from an arbitrary point set
/// \f$ P\f$ and by specialized construction methods expecting no,
/// one, two or three points as arguments. The latter methods can be
/// useful for reconstructing \f$ mc(P)\f$ from a given support set
/// \f$ S\f$ of \f$ P\f$.
/// @{

/*!

initializes `min_circle` to \f$ mc(P)\f$ with \f$ P\f$ being the set of points
in the range [`first`,`last`). If `randomize` is
`true`, a random permutation of \f$ P\f$ is computed in
advance, using the random numbers generator `random`.
Usually, this will not be necessary, however, the algorithm's
efficiency depends on the order in which the points are
processed, and a bad order might lead to extremely poor
performance (see example below).
\tparam InputIterator must be a model of `InputIterator` with `Point` as value type.
*/
template < class InputIterator >
Min_circle_2( InputIterator first,
InputIterator last,
bool randomize,
Random& random = CGAL::get_default_random(),
const Traits& traits = Traits() );

/*!

initializes `min_circle` to
\f$ mc(\emptyset)\f$, the empty set.
\post `min_circle.is_empty()` = `true`.
*/
Min_circle_2( const Traits& traits = Traits());

/*!

initializes `min_circle` to \f$ mc(\{p\})\f$, the set \f$ \{p\}\f$.
\post `min_circle.is_degenerate()` = `true`.
*/
Min_circle_2( const Point& p,
const Traits& traits = Traits());

/*!

initializes `min_circle` to \f$ mc(\{p1,p2\})\f$, the circle with diameter
equal to the segment connecting \f$ p1\f$ and \f$ p2\f$.
*/
Min_circle_2( const Point& p1,
const Point& p2,
const Traits& traits = Traits());

/*!

initializes `min_circle` to \f$ mc(\{p1,p2,p3\})\f$.
*/
Min_circle_2( const Point& p1,
const Point& p2,
const Point& p3,
const Traits& traits = Traits());

/// @}

/// \name Access Functions
/// @{

/*!

returns the number of points of `min_circle`, i.e.\ \f$ |P|\f$.
*/
int number_of_points( ) const;

/*!

returns the number of support points of `min_circle`, i.e.\ \f$ |S|\f$.
*/
int number_of_support_points( ) const;

/*!

returns an iterator referring to the first point of `min_circle`.
*/
Point_iterator points_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Point_iterator points_end() const;

/*!

returns an iterator referring to the first support point of `min_circle`.
*/
Support_point_iterator support_points_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Support_point_iterator support_points_end() const;

/*!

returns the `i`-th support point of `min_circle`. Between two
modifying operations (see below) any call to
`min_circle.support_point(i)` with the same `i` returns
the same point.
\pre \f$ 0 \leq i< \f$ `min_circle.number_of_support_points()`.
*/
const Point& support_point( int i) const;

/*!

returns the current circle of `min_circle`.
*/
const Circle& circle( ) const;

/// @}

/// \name Predicates
/// By definition, an empty `Min_circle_2` has no boundary and no
/// bounded side, i.e.\ its unbounded side equals the whole space \f$
/// \E^2\f$.
/// @{

/*!

returns `CGAL::ON_BOUNDED_SIDE`,
`CGAL::ON_BOUNDARY`, or
`CGAL::ON_UNBOUNDED_SIDE` iff `p` lies properly
inside, on the boundary of, or properly outside of `min_circle`, resp.
*/
CGAL::Bounded_side
bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies properly inside `min_circle`.
*/
bool has_on_bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies on the boundary
of `min_circle`.
*/
bool has_on_boundary( const Point& p) const;

/*!

returns `true`, iff `p` lies properly outside of `min_circle`.
*/
bool has_on_unbounded_side( const Point& p) const;

/*!

returns `true`, iff `min_circle` is empty (this implies
degeneracy).
*/
bool is_empty( ) const;

/*!

returns `true`, iff `min_circle` is degenerate,
i.e.\ if `min_circle` is empty or equal to a single point, equivalently
if the number of support points is less than 2.
*/
bool is_degenerate( ) const;

/// @}

/// \name Modifiers
/// New points can be added to an existing `min_circle`, allowing to
/// build \f$ mc(P)\f$ incrementally, e.g. if \f$ P\f$ is not known in
/// advance. Compared to the direct creation of \f$ mc(P)\f$, this is
/// not much slower, because the construction method is incremental
/// itself.
/// @{

/*!

inserts `p` into `min_circle` and recomputes the smallest
enclosing circle.
*/
void insert( const Point& p);

/*!

inserts the points in the range [`first`,`last`)
into `min_circle` and recomputes the smallest enclosing circle by
calling `insert(p)` for each point `p` in
[`first`,`last`).
\tparam InputIterator must be a model of `InputIterator` with `Point` as value type.
*/
template < class InputIterator >
void insert( InputIterator first,
InputIterator last );

/*!

deletes all points in `min_circle` and sets `min_circle` to the empty set.
\post `min_circle.is_empty()` = `true`.
*/
void clear( );

/// @}

/// \name Validity Check
/// An object `min_circle` is valid, iff <UL> <LI>`min_circle`
/// contains all points of its defining set \f$ P\f$, <LI>`min_circle`
/// is the smallest circle spanned by its support set \f$ S\f$, and
/// <LI>\f$ S\f$ is minimal, i.e.\ no support point is redundant. </UL>
/// @{

/*!

returns `true`, iff `min_circle` is valid. If `verbose`
is `true`, some messages concerning the performed checks
are written to standard error stream. The second parameter
`level` is not used, we provide it only for consistency
with interfaces of other classes.
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

}; /* end Min_circle_2 */

/*!

writes `min_circle` to output stream `os`.
An overload of `operator<<` must be defined for `Point` (and for `Circle`, if pretty printing is used).
\relates Min_circle_2
*/
std::ostream&
operator << ( std::ostream& os,
const Min_circle_2<Traits>& min_circle);

/*!

reads `min_circle` from input stream `is`.
An overload of `operator>>` must be defined for `Point`.
\relates Min_circle_2
*/
std::istream&
operator >> ( std::istream& is,
Min_circle_2<Traits> min_circle&);


} /* end namespace CGAL */

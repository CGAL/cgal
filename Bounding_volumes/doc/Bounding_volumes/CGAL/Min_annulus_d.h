
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

An object of the class `Min_annulus_d` is the unique annulus (region between
two concentric spheres with radii \f$ r\f$ and \f$ R\f$, \f$ r \leq R\f$) enclosing a
finite set of points in \f$ d\f$-dimensional Euclidean space \f$ \E^d\f$, where the
difference \f$ R^2-r^2\f$ is minimal. For a point set \f$ P\f$ we denote by \f$ ma(P)\f$
the smallest annulus that contains all points of \f$ P\f$. Note that \f$ ma(P)\f$
can be degenerate,
i.e.\ \f$ ma(P)=\emptyset\f$ if
\f$ P=\emptyset\f$ and \f$ ma(P)=\{p\}\f$ if
\f$ P=\{p\}\f$.

An inclusion-minimal subset \f$ S\f$ of \f$ P\f$ with \f$ ma(S)=ma(P)\f$ is called a
<I>support set</I>,
the points in \f$ S\f$ are the <I>support points</I>. A support set has size at
most \f$ d+2\f$, and all its points lie on the boundary of \f$ ma(P)\f$. In general,
the support set is not necessarily unique.

The underlying algorithm can cope with all kinds of input, e.g. \f$ P\f$ may be
empty or points may occur more than once. The algorithm computes a support
set \f$ S\f$ which remains fixed until the next set, insert, or clear operation.

\tparam Traits must be a model for `MinSphereAnnulusDTraits`.

We provide the models `Min_sphere_annulus_d_traits_2`,
`Min_sphere_annulus_d_traits_3`, and `Min_sphere_annulus_d_traits_d` using the
two-, three-, and \f$ d\f$-dimensional \cgal kernel, respectively.

\sa `CGAL::Min_sphere_d<Traits>`
\sa `CGAL::Min_sphere_annulus_d_traits_2<K,ET,NT>`
\sa `CGAL::Min_sphere_annulus_d_traits_3<K,ET,NT>`
\sa `CGAL::Min_sphere_annulus_d_traits_d<K,ET,NT>`
\sa `MinSphereAnnulusDTraits`

\cgalHeading{Implementation}

The problem of finding the smallest enclosing annulus of a finite point set
can be formulated as an optimization problem with linear constraints and a
linear objective
function.
The solution is obtained using our exact
solver for linear and quadratic programs \cgalCite{gs-eegqp-00}.

The creation time is almost always linear in the number of points. Access
functions and predicates take constant time, inserting a point takes almost
always linear time. The clear operation and the check for validity each
take linear time.

*/
template< typename Traits >
class Min_annulus_d {
public:

/// \name Types
/// @{

/*!
typedef to `Traits::Point_d`.
Point type used to represent the input points.
*/
typedef unspecified_type Point;

/*!
typedef to `Traits::FT`.
Number type used to return the squared radii of the smallest
enclosing annulus.
*/
typedef unspecified_type FT;

/*!
typedef to `Traits::ET`.
Number type used to do the exact computations in the underlying
solver for quadratic programs (cf. <B>Implementation</B>).
*/
typedef unspecified_type ET;

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `Point`. Used to access the points
of the smallest enclosing annulus.
*/
typedef unspecified_type Point_iterator;

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `Point`. Used to access the support points
of the smallest enclosing annulus.
*/
typedef unspecified_type Support_point_iterator;

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `Point`. Used to access the inner support points
of the smallest enclosing annulus.
*/
typedef unspecified_type Inner_support_point_iterator;

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `Point`. Used to access the outer support points
of the smallest enclosing annulus.
*/
typedef unspecified_type Outer_support_point_iterator;

/*!

non-mutable model of the \stl concept <I>RandomAccessIterator</I>
with value type `ET`. Used to access the coordinates of
the center of the smallest enclosing annulus.
*/
typedef unspecified_type Coordinate_iterator;

/// @}

/// \name Creation
/// @{

/*!

initializes `min_annulus` to \f$ ma(\emptyset)\f$.
*/
Min_annulus_d( const Traits& traits = Traits(),
int verbose = 0,
std::ostream& stream = std::cout);

/*!

initializes `min_annulus` to \f$ ma(P)\f$ with \f$ P\f$ being the set of points
in the range [`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
\pre All points have the same dimension.
*/
template < class InputIterator >
Min_annulus_d( InputIterator first,
InputIterator last,
const Traits& traits = Traits(),
int verbose = 0,
std::ostream& stream = std::cout);

/// @}

/// \name Access Functions
/// @{

/*!

returns the dimension of the points in \f$ P\f$.
If `min_annulus` is empty, the ambient dimension is \f$ -1\f$.
*/
int ambient_dimension( ) const;

/*!

returns the number of points of `min_annulus`, i.e.\ \f$ |P|\f$.
*/
int number_of_points( ) const;

/*!

returns the number of support points of `min_annulus`, i.e.\ \f$ |S|\f$.
*/
int number_of_support_points( ) const;

/*!

returns the number of support points of `min_annulus`
which lie on the inner sphere.
*/
int number_of_inner_support_points( ) const;

/*!

returns the number of support points of `min_annulus`
which lie on the outer sphere.
*/
int number_of_outer_support_points( ) const;

/*!

returns an iterator referring to the first point of `min_annulus`.
*/
Point_iterator points_begin( ) const;

/*!

returns the corresponding past-the-end iterator.
*/
Point_iterator points_end( ) const;

/*!

returns an iterator referring to the first support point of `min_annulus`.
\pre \f$ ma(P)\f$ is not degenerate, i.e., `number_of_support_points()` is at least one.
*/
Support_point_iterator support_points_begin( ) const;

/*!

returns the corresponding past-the-end iterator.
\pre \f$ ma(P)\f$ is not degenerate, i.e., `number_of_support_points()` is at least one.
*/
Support_point_iterator support_points_end( ) const;

/*!

returns an iterator referring to the first inner support point
of `min_annulus`.
*/
Inner_support_point_iterator
inner_support_points_begin( ) const;

/*!

returns the corresponding past-the-end iterator.
*/
Inner_support_point_iterator
inner_support_points_end( ) const;

/*!

returns an iterator referring to the first outer support point
of `min_annulus`.
*/
Outer_support_point_iterator
outer_support_points_begin( ) const;

/*!

returns the corresponding past-the-end iterator.
*/
Outer_support_point_iterator
outer_support_points_end( ) const;

/*!

returns the center of `min_annulus`.
An implicit conversion from `ET` to `RT` must be available.
\pre `min_annulus` is not empty.
*/
Point center( ) const;

/*!

returns the squared inner radius of `min_annulus`.
An implicit conversion from `ET` to `RT` must be available.
\pre `min_annulus` is not empty.
*/
FT squared_inner_radius( ) const;

/*!

returns the squared outer radius of `min_annulus`.
An implicit conversion from `ET` to `RT` must be available.
\pre `min_annulus` is not empty.
*/
FT squared_outer_radius( ) const;

/*!

returns an iterator referring to the first coordinate
of the center of `min_annulus`.

\note
The coordinates have a rational
representation, i.e.\ the first \f$ d\f$ elements of the iterator
range are the numerators and the \f$ (d\!+\!1)\f$-st element is the
common denominator.
*/
Coordinate_iterator
center_coordinates_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Coordinate_iterator
center_coordinates_end() const;

/*!

returns the numerator of the squared inner radius of `min_annulus`.
*/
ET squared_inner_radius_numerator( ) const;

/*!

returns the numerator of the squared outer radius of `min_annulus`.
*/
ET squared_outer_radius_numerator( ) const;

/*!

returns the denominator of the squared radii of `min_annulus`.
*/
ET squared_radii_denominator( ) const;

/// @}

/// \name Predicates
/// The bounded area of the smallest enclosing annulus lies between
/// the inner and the outer sphere. The boundary is the union of both
/// spheres. By definition, an empty annulus has no boundary and no
/// bounded side, i.e.\ its unbounded side equals the whole space \f$
/// \E^d\f$.
/// @{

/*!

returns `CGAL::ON_BOUNDED_SIDE`,
`CGAL::ON_BOUNDARY`, or
`CGAL::ON_UNBOUNDED_SIDE` iff `p` lies
properly inside, on the boundary, or properly outside of
`min_annulus`, resp.
\pre The dimension of `p` equals `min_annulus.ambient_dimension()` if `min_annulus` is not empty.
*/
CGAL::Bounded_side
bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies properly inside `min_annulus`.
\pre The dimension of `p` equals `min_annulus.ambient_dimension()` if `min_annulus` is not empty.
*/
bool has_on_bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies on the boundary of `min_annulus`.
\pre The dimension of `p` equals `min_annulus.ambient_dimension()` if `min_annulus` is not empty.
*/
bool has_on_boundary( const Point& p) const;

/*!

returns `true`, iff `p` lies properly outside of `min_annulus`.
\pre The dimension of `p` equals `min_annulus.ambient_dimension()` if `min_annulus` is not empty.
*/
bool has_on_unbounded_side( const Point& p) const;

/*!

returns `true`, iff `min_annulus` is empty (this implies degeneracy).
*/
bool is_empty( ) const;

/*!

returns `true`, iff `min_annulus` is degenerate, i.e.\ if `min_annulus` is empty or equal to a single point.
*/
bool is_degenerate( ) const;

/// @}

/// \name Modifiers
/// @{

/*!

resets `min_annulus` to \f$ ma(\emptyset)\f$.
*/
void clear( );

/*!

sets `min_annulus` to \f$ ma(P)\f$, where \f$ P\f$ is the set of points in
the range [`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
\pre All points have the same dimension.
*/
template < class InputIterator >
void set( InputIterator first,
InputIterator last );

/*!

inserts `p` into `min_annulus`.
\pre The dimension of `p` equals `min_annulus.ambient_dimension()` if `min_annulus` is not empty.
*/
void insert( const Point& p);

/*!

inserts the points in the range [`first`,`last`) into
`min_annulus` and recomputes the smallest enclosing annulus.
@tparam InputIterator is a model of `InputIterator` with `Point` as value type.
\pre All points have the same dimension. If `min_annulus` is not empty, this dimension must be equal to `min_annulus.ambient_dimension()`.
*/
template < class InputIterator >
void insert( InputIterator first,
InputIterator last );

/// @}

/// \name Validity Check
/// An object `min_annulus` is valid, iff <UL> <LI>`min_annulus`
/// contains all points of its defining set \f$ P\f$,
/// <LI>`min_annulus` is the smallest annulus containing its support
/// set \f$ S\f$, and <LI>\f$ S\f$ is minimal, i.e.\ no support point
/// is redundant. </UL> <I>Note:</I> In this release only the first
/// item is considered by the validity check.
/// @{

/*!

returns `true`, iff `min_annulus` is valid. If `verbose` is
`true`, some messages concerning the performed checks are
written to standard error stream. The second parameter
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

}; /* end Min_annulus_d */

/*!

writes `min_annulus` to output stream `os`.
An overload of `operator<<` must be defined for `Point`.
\relates Min_annulus_d
*/
std::ostream&
operator << ( std::ostream& os,
const Min_annulus_d<Traits>& min_annulus);

/*!

reads `min_annulus` from input stream `is`.
An overload of `operator>>` must be defined for `Point`.
\relates Min_annulus_d
*/
std::istream&
operator >> ( std::istream& is,
Min_annulus_d<Traits> min_annulus&);


} /* end namespace CGAL */

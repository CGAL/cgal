
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

An object of the class `Min_sphere_d` is the unique sphere of
smallest volume enclosing a finite (multi)set of points in \f$ d\f$-dimensional
Euclidean space \f$ \E^d\f$. For a set \f$ P\f$ we denote by \f$ ms(P)\f$ the
smallest sphere that contains all points of \f$ P\f$. \f$ ms(P)\f$ can
be degenerate, i.e.\ \f$ ms(P)=\emptyset\f$
if \f$ P=\emptyset\f$ and \f$ ms(P)=\{p\}\f$ if
\f$ P=\{p\}\f$.

An inclusion-minimal subset \f$ S\f$ of \f$ P\f$ with \f$ ms(S)=ms(P)\f$ is called a
<I>support set</I>, the points in \f$ S\f$ are the <I>support points</I>.
A support set has size at most \f$ d+1\f$, and all its points lie on the
boundary of \f$ ms(P)\f$. In general, neither the support set nor its size
are unique.

The algorithm
computes a support set \f$ S\f$ which remains fixed until the next insert
or clear operation.

\note  This class is (almost) obsolete. The class
`CGAL::Min_sphere_of_spheres_d<Traits>` solves a more general problem
and is faster than `Min_sphere_d` even if used only for points
as input. Most importantly, `CGAL::Min_sphere_of_spheres_d<Traits>` has
a specialized implementation for floating-point arithmetic which
ensures correct results in a large number of cases (including
highly degenerate ones). In contrast, `Min_sphere_d` is not reliable
under floating-point computations. The only advantage of
`Min_sphere_d` over `CGAL::Min_sphere_of_spheres_d<Traits>` is that the
former can deal with points in homogeneous coordinates, in which
case the algorithm is division-free. Thus, `Min_sphere_d` might still
be an option in case your input number type cannot (efficiently)
divide.


\tparam Traits must be a model of the concept
`MinSphereAnnulusDTraits` as its template argument.

We provide the models `CGAL::Min_sphere_annulus_d_traits_2`,
`CGAL::Min_sphere_annulus_d_traits_3` and
`CGAL::Min_sphere_annulus_d_traits_d`
for two-, three-, and \f$ d\f$-dimensional points respectively.

\sa `CGAL::Min_sphere_annulus_d_traits_2<K,ET,NT>`
\sa `CGAL::Min_sphere_annulus_d_traits_3<K,ET,NT>`
\sa `CGAL::Min_sphere_annulus_d_traits_d<K,ET,NT>`
\sa `MinSphereAnnulusDTraits`
\sa `CGAL::Min_circle_2<Traits>`
\sa `CGAL::Min_sphere_of_spheres_d<Traits>`
\sa `CGAL::Min_annulus_d<Traits>`

\cgalHeading{Implementation}

We implement the algorithm of Welzl with move-to-front
heuristic \cgalCite{w-sedbe-91a} for small point sets, combined with a new
efficient method for large sets, which is particularly tuned for
moderately large dimension (\f$ d \leq 20\f$) \cgalCite{cgal:g-frseb-99}.
The creation time is almost
always linear in the number of points. Access functions and predicates
take constant time, inserting a point might take up to linear time,
but substantially less than computing the new smallest enclosing
sphere from scratch. The clear operation and the check for validity
each take linear time.

\cgalHeading{Example}

\cgalExample{Min_sphere_d/min_sphere_homogeneous_3.cpp}

*/
template< typename Traits >
class Min_sphere_d {
public:

/// \name Types
/// @{

/*!

*/
typedef unspecified_type Traits;

/*!
typedef to `Traits::FT`.
*/
typedef unspecified_type FT;

/*!
typedef to `Traits::Point`.
*/
typedef unspecified_type Point;

/*!
non-mutable model of the STL
concept `BidirectionalIterator` with value type `Point`. Used
to access the points used to build the smallest enclosing sphere.
*/
typedef unspecified_type Point_iterator;

/*!
non-mutable model of the STL
concept `BidirectionalIterator` with value type `Point`. Used
to access the support points defining the smallest enclosing sphere.
*/
typedef unspecified_type Support_point_iterator;

/// @}

/// \name Creation
/// @{

/*!

creates a variable of type `Min_sphere_d` and
initializes it to \f$ ms(\emptyset)\f$.
If the traits parameter is not supplied, the class `Traits`
must provide a default constructor.
*/
Min_sphere_d (const Traits& traits = Traits());

/*!

creates a variable `min_sphere` of type `Min_sphere_d`.
It is initialized to \f$ ms(P)\f$ with \f$ P\f$ being the set of points
in the range [`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
If the traits parameter is not supplied, the class `Traits` must provide a default constructor.
\pre All points have the same dimension.
*/
template < class InputIterator >
Min_sphere_d( InputIterator first,
InputIterator last,
const Traits& traits = Traits());

/*!

returns the number of points of `min_sphere`, i.e.\ \f$ |P|\f$.
*/
int number_of_points( ) const;

/*!

returns the number of support points of `min_sphere`, i.e.\ \f$ |S|\f$.
*/
int number_of_support_points( ) const;

/*!

returns an iterator referring to the first point of `min_sphere`.
*/
Point_iterator points_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Point_iterator points_end() const;

/*!

returns an iterator referring to the first support point of `min_sphere`.
*/
Support_point_iterator support_points_begin() const;

/*!

returns the corresponding past-the-end iterator.
*/
Support_point_iterator support_points_end() const;

/*!

returns the dimension of the points in \f$ P\f$. If `min_sphere` is empty, the ambient dimension is \f$ -1\f$.
*/
int ambient_dimension() const;

/*!

returns the center of `min_sphere`.
\pre `min_sphere` is not empty.
*/
const Point& center( ) const;

/*!

returns the squared radius of `min_sphere`.
\pre `min_sphere` is not empty.
*/
FT squared_radius( ) const;

/// @}

/// \name Predicates
/// By definition, an empty `Min_sphere_d` has no boundary and no
/// bounded side, i.e.\ its unbounded side equals the whole space \f$
/// \E^d\f$.
/// @{

/*!

returns `CGAL::ON_BOUNDED_SIDE`, `CGAL::ON_BOUNDARY`, or
`CGAL::ON_UNBOUNDED_SIDE` iff `p` lies properly inside,
on the boundary, or properly outside of `min_sphere`, resp.
\pre If `min_sphere` is not empty, the dimension of \f$ p\f$ equals `ambient_dimension()`.
*/
Bounded_side
bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies properly inside `min_sphere`.
\pre If `min_sphere` is not empty, the dimension of \f$ p\f$ equals `ambient_dimension()`.
*/
bool has_on_bounded_side( const Point& p) const;

/*!

returns `true`, iff `p` lies on the boundary
of `min_sphere`.
\pre if `min_sphere` is not empty, the dimension of \f$ p\f$ equals `ambient_dimension()`.
*/
bool has_on_boundary( const Point& p) const;

/*!

returns `true`, iff `p` lies properly outside of `min_sphere`.
\pre If `min_sphere` is not empty, the dimension of \f$ p\f$ equals `ambient_dimension()`.
*/
bool has_on_unbounded_side( const Point& p) const;

/*!

returns `true`, iff `min_sphere` is empty (this implies
degeneracy).
*/
bool is_empty( ) const;

/*!

returns `true`, iff `min_sphere` is degenerate, i.e.\ if
`min_sphere` is empty or equal to a single point, equivalently if
the number of support points is less than 2.
*/
bool is_degenerate( ) const;

/// @}

/// \name Modifiers
/// @{

/*!

resets `min_sphere` to \f$ ms(\emptyset)\f$.
*/
void clear ();

/*!

sets `min_sphere` to the \f$ ms(P)\f$, where \f$ P\f$ is the set of points
in the range [`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
\pre All points have the same dimension.
*/
template < class InputIterator >
void set( InputIterator first,
InputIterator last );

/*!

inserts `p` into `min_sphere`. If `p` lies inside the
current sphere, this is a constant-time operation, otherwise
it might take longer, but usually substantially less than
recomputing the smallest enclosing sphere from scratch.
\pre The dimension of `p` equals `ambient_dimension()` if `min_sphere` is not empty.
*/
void insert( const Point& p);

/*!

inserts the points in the range [`first`,`last`)
into `min_sphere` and recomputes the smallest enclosing sphere, by
calling `insert` for all points in the range.
\tparam InputIterator is a model of `InputIterator` with `Point` as value type.
\pre All points have the same dimension. If `min_sphere` is not empty, this dimension must be equal to `ambient_dimension()`.
*/
template < class InputIterator >
void insert( InputIterator first,
InputIterator last );

/// @}

/// \name Validity Check
/// An object `min_sphere` is valid, iff <UL> <LI>`min_sphere`
/// contains all points of its defining set \f$ P\f$, <LI>`min_sphere`
/// is the smallest sphere containing its support set \f$ S\f$, and
/// <LI>\f$ S\f$ is minimal, i.e.\ no support point is redundant. </UL>
///
/// \note Under inexact arithmetic, the result of the
/// validation is not realiable, because the checker itself can suffer
/// from numerical problems.
/// @{

/*!

returns `true`, iff `min_sphere` is valid. If `verbose`
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

}; /* end Min_sphere_d */


/*!

writes `min_sphere` to output stream `os`.
An overload of `operator<<` must be defined for `Point`.
\relates Min_sphere_d
*/
std::ostream& operator << ( std::ostream& os,
const Min_sphere_d<Traits>&
min_sphere);

/*!

reads `min_sphere` from input stream `is`.
An overload of `operator>>` must be defined for `Point`
\relates Min_sphere_d
*/
std::istream& operator >> ( std::istream& is,
Min_sphere_d<Traits> min_sphere&);

} /* end namespace CGAL */

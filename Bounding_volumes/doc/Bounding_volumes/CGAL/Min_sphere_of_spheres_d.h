
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

An object of the class `Min_sphere_of_spheres_d` is a data structure that represents
the unique sphere of smallest volume enclosing a finite set of spheres
in \f$ d\f$-dimensional Euclidean space \f$ \E^d\f$. For a set \f$ S\f$ of spheres
we denote by \f$ ms(S)\f$ the smallest sphere that contains all spheres of
\f$ S\f$; we call \f$ ms(S)\f$ the <I>minsphere</I> of \f$ S\f$. \f$ ms(S)\f$ can be
degenerate, i.e., \f$ ms(S)=\emptyset\f$,
if \f$ S=\emptyset\f$ and \f$ ms(S)=\{s\}\f$,
if \f$ S=\{s\}\f$. Any sphere in \f$ S\f$ may be degenerate, too, i.e., any
sphere from \f$ S\f$ may be a point. Also, \f$ S\f$ may contain several
copies of the same sphere.

An inclusion-minimal subset \f$ R\f$ of \f$ S\f$ with \f$ ms(R)=ms(S)\f$ is called a
<I>support set</I> for \f$ ms(S)\f$; the spheres in \f$ R\f$ are the
<I>support spheres</I>. A support set has size at most \f$ d+1\f$, and all
its spheres lie on the boundary of \f$ ms(S)\f$. (A sphere \f$ s'\f$ is said to
<I>lie on the boundary</I> of a sphere \f$ s\f$, if \f$ s'\f$ is contained in \f$ s\f$
and if their boundaries intersect.) In general, the support set is
not unique.

The algorithm computes the center and the radius of \f$ ms(S)\f$, and finds
a support set \f$ R\f$ (which remains fixed until the next `insert()`,
`clear()` or `set()` operation). We also provide a
specialization of the algorithm for the case when the center
coordinates and radii of the input spheres are floating-point numbers.
This specialized algorithm uses floating-point arithmetic only, is
very fast and especially tuned for stability and robustness. Still,
it's output may be incorrect in some (rare) cases; termination is
guaranteed.

When default constructed, an instance of type
`Min_sphere_of_spheres_d<Traits>` represents the set
\f$ S=\emptyset\f$, together with its minsphere \f$ ms(S)=\emptyset\f$. You can
add spheres to the set \f$ S\f$ by calling `insert()`. Querying the
minsphere is done by calling the routines `is_empty()`,
`radius()` and `center_cartesian_begin()`, among others.

In general, the radius and the Euclidean center coordinates of \f$ ms(S)\f$
need not be rational. Consequently, the algorithm computing the exact
minsphere will have to deal with algebraic numbers. Fortunately, both
the radius and the coordinates of the minsphere are numbers of the
form \f$ a_i+b_i\sqrt{t}\f$, where \f$ a_i,b_i,t\in \Q\f$ and where \f$ t\ge 0\f$ is
the same for all coordinates and the radius. Thus, the exact
minsphere can be described by the number \f$ t\f$, which is called the
sphere's <I>discriminant</I>, and by \f$ d+1\f$ pairs \f$ (a_i,b_i)\in\Q^2\f$
(one for the radius and \f$ d\f$ for the center coordinates).

\note This class (almost) replaces
`CGAL::Min_sphere_d<Traits>`, which solves the less general
problem of finding the smallest enclosing ball of a set of
<I>points</I>. `Min_sphere_of_spheres_d` is faster than
`CGAL::Min_sphere_d<Traits>`, and in contrast to the latter
provides a specialized implementation for floating-point arithmetic
which ensures correct results in a large number of cases (including
highly degenerate ones). The only advantage of
`CGAL::Min_sphere_d<Traits>` over `Min_sphere_of_spheres_d` is that the former
can deal with points in homogeneous coordinates, in which case the
algorithm is division-free. Thus, `CGAL::Min_sphere_d<Traits>`
might still be an option in case your input number type cannot
(efficiently) divide.

\tparam Traits must be a model of the concept
`MinSphereOfSpheresTraits` as its template argument.

\sa `CGAL::Min_sphere_d<Traits>`
\sa `CGAL::Min_circle_2<Traits>`

\cgalHeading{Implementation}

We implement two algorithms, the LP-algorithm and a
heuristic \cgalCite{msw-sblp-92}. As described in the documentation of
concept `MinSphereOfSpheresTraits`, each has its advantages and
disadvantages: Our implementation of the LP-algorithm has maximal
expected running time \cgalBigO{2^d n}, while the heuristic comes without
any complexity guarantee. In particular, the LP-algorithm runs in
linear time for fixed dimension \f$ d\f$. (These running times hold for the
arithmetic model, so they count the number of operations on
the number type `FT`.)

On the other hand, the LP-algorithm is, for inexact number types
`FT`, much worse at handling degeneracies and should therefore not
be used in such a case. (For exact number types
`FT`, both methods handle all kinds of degeneracies.)

Currently, we require `Traits::FT` to be either an exact number
type or `double` or `float`; other inexact number types are
not supported at this time. Also, the current implementation only
handles spheres with %Cartesian coordinates; homogeneous representation
is not supported yet.

\cgalHeading{Example}

\cgalExample{Min_sphere_of_spheres_d/min_sphere_of_spheres_d_d.cpp}

*/
template< typename Traits >
class Min_sphere_of_spheres_d {
public:

/// \name Types
/// @{

/*!
is a typedef to `Traits::Sphere`.
*/
typedef unspecified_type Sphere;

/*!
is a typedef to `Traits::FT`.
*/
typedef unspecified_type FT;

/*!
is the type of the radius and of the center
coordinates of the computed minsphere: When `FT` is an inexact
number type (`double`, for instance), then `Result` is
simply `FT`. However, when `FT` is an exact number type,
then `Result` is a typedef to a derived class of
`std::pair<FT,FT>`; an instance of this type represents the
number \f$ a+b\sqrt{t}\f$, where \f$ a\f$ is the first and \f$ b\f$ the second
element of the pair and where the number \f$ t\f$ is accessed using the
member function `discriminant()` of class
`Min_sphere_of_spheres_d<Traits>`.
*/
typedef unspecified_type Result;

/*!
is either `CGAL::LP_algorithm` or
`CGAL::Farthest_first_heuristic`. As is described in the
documentation of concept `MinSphereOfSpheresTraits`, the type
`Algorithm` reflects the method which is used to compute the
minsphere. (Normally, `Algorithm` coincides with
`Traits::Algorithm`. However, if the method
`Traits::Algorithm` should not be supported anymore in a future
release, then `Algorithm` will have another type.)
*/
typedef unspecified_type Algorithm;

/*!
non-mutable model of the STL
concept
`BidirectionalIterator`
with value type `Sphere`. Used
to access the support spheres defining the smallest enclosing sphere.
*/
typedef unspecified_type Support_iterator;

/*!
non-mutable model of the STL
concept
`BidirectionalIterator`
to access the center coordinates of the minsphere.
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

/// \name Creation
/// @{

/*!
creates a variable of type `Min_sphere_of_spheres_d` and initializes it to
\f$ ms(\emptyset)\f$. If the traits
parameter is not supplied, the class `Traits` must provide a
default constructor.
*/
Min_sphere_of_spheres_d(const Traits& traits = Traits());

/*!
creates a variable `minsphere` of type
`Min_sphere_of_spheres_d` and inserts (cf.
`insert()`) the spheres from
the range [`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Sphere` as value type.
If the traits parameter is not supplied, the class `Traits` must provide a default constructor.
*/
template < typename InputIterator >
Min_sphere_of_spheres_d( InputIterator first,
InputIterator last,
const Traits& traits = Traits());

/// @}

/// \name Access Functions
/// @{

/*!
returns
an iterator referring to the first support sphere of `minsphere`.
*/
Support_iterator support_begin() const;

/*!
returns the
corresponding past-the-end iterator.
*/
Support_iterator support_end() const;

/*!
returns the
discriminant of `minsphere`. This number is undefined when `FT` is
an inexact number type. When `FT` is exact, the center
coordinates and the radius of the minsphere are numbers of the form
\f$ a+b\sqrt{t}\f$, where \f$ t\f$ is the discriminant of the minsphere as
returned by this function. \pre `minsphere` is not empty, and `FT` is an exact number type.
*/
const FT& discriminant( ) const;

/*!
returns the radius of
`minsphere`. If `FT` is an exact number type then the radius of the
minsphere is the real number \f$ a+b\sqrt{t}\f$, where \f$ t\f$ is the
minsphere's discriminant, \f$ a\f$ is the first and \f$ b\f$ the second
component of the pair returned by `radius()`.
\pre `minsphere` is not empty.
*/
Result radius( ) const;

/*!
returns a const-iterator to the first of the
`Traits::%D` center coordinates of `minsphere`. The iterator returns
objects of type `Result`. If `FT` is an exact number type,
then a center coordinate is represented by a pair \f$ (a,b)\f$ describing
the real number \f$ a+b\sqrt{t}\f$, where \f$ t\f$ is the minsphere's
discriminant (cf. `discriminant()`). \pre `minsphere` is not empty.
*/
Cartesian_const_iterator center_cartesian_begin( )
const;

/*!
returns the corresponding past-the-end iterator, i.e.\
`center_cartesian_begin()+Traits::%D`.
\pre `minsphere` is not empty.
*/
Cartesian_const_iterator center_cartesian_end( )
const;

/// @}

/// \name Predicates
/// @{

/*!
returns `true`, iff
`minsphere` is empty, i.e.\ iff \f$ ms(S)=\emptyset\f$.
*/
bool is_empty( ) const;

/// @}

/// \name Modifiers
/// @{

/*!

resets `minsphere` to \f$ ms(\emptyset)\f$, with \f$ S:= \emptyset\f$.
*/
void clear ();

/*!
sets `minsphere` to
the \f$ ms(S)\f$, where \f$ S\f$ is the set of spheres in the range
[`first`,`last`).
\tparam InputIterator is a model of `InputIterator` with `Sphere` as value type.
*/
template < class InputIterator > void
set( InputIterator first, InputIterator last );

/*!
inserts the
sphere `s` into the set \f$ S\f$ of instance `minsphere`.
*/
void insert( const Sphere& s );

/*!
inserts the spheres in
the range [`first`,`last`) into the set \f$ S\f$ of instance
`minsphere`.
\tparam InputIterator is a model of `InputIterator` with `Sphere` as value type.
*/
template < class InputIterator > void insert(
InputIterator first, InputIterator last );

/// @}

/// \name Validity Check
/// An object `minsphere` is valid, iff <UL> <LI>`minsphere` contains
/// all spheres of its defining set \f$ S\f$, <LI>`minsphere` is the
/// smallest sphere containing its support set \f$ R\f$, and <LI>\f$
/// R\f$ is minimal, i.e., no support sphere is redundant. </UL>
/// @{

/*!
returns `true`, iff
`minsphere` is valid. When `FT` is inexact, this routine always
returns `true`.
*/
bool is_valid() const;

/// @}

/// \name Miscellaneous
/// @{

/*!

returns a const reference to the traits class object.
*/
const Traits& traits( ) const;

/// @}

}; /* end Min_sphere_of_spheres_d */
} /* end namespace CGAL */

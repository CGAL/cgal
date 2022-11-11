
/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

A model of concept `MinSphereOfSpheresTraits` must provide the
following constants, types, predicates and operations.

\cgalHasModel `CGAL::Min_sphere_of_spheres_d_traits_2<K,FT,UseSqrt,Algorithm>`
\cgalHasModel `CGAL::Min_sphere_of_spheres_d_traits_3<K,FT,UseSqrt,Algorithm>`
\cgalHasModel `CGAL::Min_sphere_of_spheres_d_traits_d<K,FT,Dim,UseSqrt,Algorithm>`

\cgalHasModel `CGAL::Min_sphere_of_points_d_traits_2<K,FT,UseSqrt,Algorithm>`
\cgalHasModel `CGAL::Min_sphere_of_points_d_traits_3<K,FT,UseSqrt,Algorithm>`
\cgalHasModel `CGAL::Min_sphere_of_points_d_traits_d<K,FT,Dim,UseSqrt,Algorithm>`
*/

class MinSphereOfSpheresTraits {
public:

/// \name Constants
/// @{

/*!
specifies the dimension of the spheres you want to
compute the minsphere of.
*/
static const int D;

/// @}

/// \name Types
/// @{

/*!
is a typedef to to some class representing a sphere.
(The package will compute the minsphere of spheres of type
`Sphere`.) The type `Sphere` must provide a copy constructor.
*/
typedef unspecified_type Sphere;

/*!
is a (exact or inexact) field number type.

\tparam FT must either be `double` or `float`, or an exact field number type. (An <I>exact</I> number type is one which evaluates arithmetic expressions involving the four basic operations and comparisions with infinite precision, that is, like in \f$ \mathbb{R}\f$.)
*/
typedef unspecified_type FT;

/*!
non-mutable model of the STL
concept `ForwardIterator` with value type `FT`. Used
to access the center coordinates of a sphere.
*/
typedef unspecified_type Cartesian_const_iterator;

/*!
must typedef to either
`CGAL::Tag_true` or `CGAL::Tag_false`. The algorithm
uses (depending on the type
`MinSphereOfSpheresTraits::Algorithm`) floating-point arithmetic
internally for some intermediate computations. The type
`Use_square_roots` affects how these calculations are done: When
`Use_square_roots` is `Tag_true`, the algorithm computing
the minsphere will perform square-root operations on `double`s
and `float`s where appropriate. On the other hand, if
`Use_square_roots` is `CGAL::Tag_false`, the algorithm will
work without doing square-roots.

<I>Note:</I> On some platforms the algorithm is much faster when
square-roots are disabled (due to lacking hardware support).
*/
typedef unspecified_type Use_square_roots;

/*!
selects the method to compute the minsphere
with. It must typedef to either `CGAL::Default_algorithm`,
`CGAL::LP_algorithm` or `CGAL::Farthest_first_heuristic`.
The recommended choice is the first, which is a synonym to the one
of the other two methods which we consider "the best in practice."
In case of `CGAL::LP_algorithm`, the minsphere will be computed
using the LP-algorithm \cgalCite{msw-sblp-92}, which in our
implementation has maximal expected running time \f$ O(2^d n)\f$ (in the
number of operations on the number type `FT`). In case of
`CGAL::Farthest_first_heuristic`, a simple heuristic will be
used instead which seems to work fine in practice, but comes without
a guarantee on the running time. For an inexact number type
`FT` we strongly recommend `CGAL::Default_algorithm`, or, if
you want, `CGAL::Farthest_first_heuristic`, since these handle
most degeneracies in a satisfying manner. Notice that this
compile-time flag is taken as a hint only. Should one of the
methods not be available anymore in a future release, then the
default algorithm will be chosen.
*/
typedef unspecified_type Algorithm;

/// @}

/// \name Access Functions
/// @{

/*!
returns the radius of sphere `s`.
\post The returned number is greater or equal to \f$ 0\f$.
*/
FT radius(const
Sphere& s);

/*!

returns an iterator referring to the first of the `D` %Cartesian
coordinates of the center of `s`.
*/
Cartesian_const_iterator center_cartesian_begin(const Sphere& s);

/// @}

}; /* end MinSphereOfSpheresTraits */

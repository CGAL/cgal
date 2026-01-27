/*!
\ingroup PkgConvexHullDConcepts
\cgalConcept

\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead.

Requirements of the traits class to be used with the
class `CGAL::Convex_hull_d`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Cartesian_d<FT,LA>}
\cgalHasModels{CGAL::Homogeneous_d<RT,LA>}
\cgalHasModels{CGAL::Convex_hull_d_traits_3<R>}
\cgalHasModelsEnd

*/
class ConvexHullTraits_d {
public:

/// \name Types
/// @{

/*!
the dD point type on which the convex hull algorithm
operates
*/
typedef unspecified_type Point_d;

/*!
a dD plane
*/
typedef unspecified_type Hyperplane_d;

/*!
a dD vector
*/
typedef unspecified_type Vector_d;

/*!
a dD ray
*/
typedef unspecified_type Ray_d;

/*!
an arithmetic ring type
*/
typedef unspecified_type RT;

/*!
Function object type that provides
`Vector_d operator()(int d, CGAL::Null_vector)`, which constructs
and returns the null vector.
*/
typedef unspecified_type Construct_vector_d;

/*!
Function object type that
provides `Hyperplane_d operator()(ForwardIterator first,
ForwardIterator last, Point_d p, CGAL::Oriented_side side)`, which
constructs and returns a hyperplane passing through the points in
`tuple[first,last)` and oriented such that `p` is on the side
`side` of the returned hyperplane. When
`side==ON_ORIENTED_BOUNDARY` then any hyperplane containing the
tuple is returned.
*/
typedef unspecified_type Construct_hyperplane_d;

/*!
Function object type that provides
`Point_d operator()(Vector_d v)`, which constructs and
returns the point defined by \f$ 0+v\f$.
*/
typedef unspecified_type Vector_to_point_d;

/*!
Function object type that provides
`Vector_d operator()(Point_d v)`, which constructs and returns the
vector defined by \f$ p-0\f$.
*/
typedef unspecified_type Point_to_vector_d;

/*!
Function object type that provides
`Orientation operator()(ForwardIterator first,
ForwardIterator last)`, which determines the orientation of the
points `tuple[first,last)`.
*/
typedef unspecified_type Orientation_d;

/*!
Function object type that provides
`Vector_d operator()(Hyperplane_d h)`, which constructs and
returns a vector orthogonal to `h` and pointing from the boundary
into its positive halfspace.
*/
typedef unspecified_type Orthogonal_vector_d;

/*!
Predicate object type that provides
`Oriented_side operator()(Hyperplane_d h, Point_d p)`, which
determines the oriented side of `p` with respect to `h`.
*/
typedef unspecified_type Oriented_side_d;

/*!
Predicate object type that
provides `bool operator()(Hyperplane_d h, Point_d p)`, which
return true iff `p` lies in the positive halfspace determined by
`h`.
*/
typedef unspecified_type Has_on_positive_side_d;

/*!
Predicate object type that provides
`bool operator()(ForwardIterator first, ForwardIterator last)`, which
determines if the points `tuple[first,last)` are affinely independent.
*/
typedef unspecified_type Affinely_independent_d;

/*!
Predicate object type that
provides `bool operator()(ForwardIterator first, ForwardIterator
last, Point_d p)`, which determines if `p` is contained in
the closed simplex defined by the points in `tuple[first,last)`.
*/
typedef unspecified_type Contained_in_simplex_d;

/*!
Predicate object type that
provides `bool operator()(ForwardIterator first, ForwardIterator
last, Point_d p)`, which determines if `p` is contained in
the affine hull of the points in `tuple[first,last)`.
*/
typedef unspecified_type Contained_in_affined_hull_d;

/*!
Predicate object type that provides
`Object operator()(Ray_d r, Hyperplane_d h)`, which determines if
`r` and `h` intersect and returns the corresponding
polymorphic object.
*/
typedef unspecified_type Intersect_d;

/// @}

/// A default constructor and copy constructor is required.
ConvexHullTraits_d();

/*! \name Operations
For each of the above function and predicate object types,
`Func_obj_type`, a function must exist with the name
`func_obj_type_object` that creates an instance of the function or
predicate object type. For example:
*/
/// @{

/*!

*/
Construct_vector_d construct_vector_d_object();

/// @}

}; /* end ConvexHullTraits_d */

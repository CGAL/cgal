/*!
\ingroup PkgConvexHullDConcepts
\cgalConcept

\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead.

Requirements of the first traits class to be used with the
class `CGAL::Delaunay_d`.

\cgalHasModel `CGAL::Cartesian_d<FT,LA>`
\cgalHasModel `CGAL::Homogeneous_d<RT,LA>`

*/

class DelaunayTraits_d {
public:

/// \name Types
/// @{

/*!
the dD point type on which the Delaunay
algorithm operates
*/
typedef unspecified_type Point_d;

/*!
a dD sphere
*/
typedef unspecified_type Sphere_d;

/*!
an arithmetic field type
*/
typedef unspecified_type FT;

/*!
Predicate object type that provides
`Point_d operator()(Sphere_d s, int i)`, which returns the \f$ i\f$th
point defining sphere `s`.
*/
typedef unspecified_type Point_of_sphere_d;

/*!
Predicate object type that provides
`Sphere_d operator()(int d, ForwardIterator first, ForwardIterator
last)`, which returns a dD sphere through the points in
`tuple[first,last)`.
*/
typedef unspecified_type Construct_sphere_d;

/*!
Predicate object type that
provides `bool operator()(ForwardIterator first, ForwardIterator
last, Point_d p)`, which determines if `p` is contained in
the closed simplex defined by the points in `tuple[first,last)`.
*/
typedef unspecified_type Contained_in_simplex_d;

/*!
Predicate object type that provides
`FT operator()(Point_d p,Point_d q)`, which determines the
squared distance from `p` to `q`.
*/
typedef unspecified_type Squared_distance_d;

/*!
Predicate object type that
provides `bool operator()(ForwardIterator first, ForwardIterator
last)`, which determines if the points in `tuple[first,last)` are
affinely independent.
*/
typedef unspecified_type Affinely_independent_d;

/// @}

/// A default constructor and copy constructor is required.
DelaunayTraits_d();

/*! \name Operations
For each of the above function and predicate object types,
`Func_obj_type`, a function must exist with the name
`func_obj_type_object` that creates an instance of the function or
predicate object type. For example:
*/
/// @{

/*!

*/
Construct_sphere_d construct_sphere_d_object();

/// @}

}; /* end DelaunayTraits_d */

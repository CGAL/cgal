
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

In addition to the requirements of the concept `TriangulationTraits_2`
the concept
`DelaunayTriangulationTraits_2` requires a predicate to check the empty
circle property. The corresponding predicate type is called type
`Side_of_oriented_circle_2`.

The additional types `Line_2`,
`Ray_2` and the constructor objects
`Construct_ray_2`,
`Construct_circumcenter_2`, `Construct_bisector_2`,
`Construct_midpoint`
are used to build the dual Voronoi diagram
and are required only if the dual functions are called.
The additional predicate type `Compare_distance_2` is
required if the  method `nearest_vertex()` is used.

\cgalRefines `TriangulationTraits_2`


\cgalHasModel \cgal kernels
\cgalHasModel `CGAL::Projection_traits_xy_3<K>` (not for dual Voronoi functions)
\cgalHasModel `CGAL::Projection_traits_yz_3<K>` (not for dual Voronoi functions)
\cgalHasModel `CGAL::Projection_traits_xz_3<K>` (not for dual Voronoi functions)

\sa `TriangulationTraits_2`
*/

class DelaunayTriangulationTraits_2 {
public:

/// \name Types
/// @{

/*!
The line type. This type is required only if
some dual functions are called.
*/
typedef unspecified_type Line_2;

/*!
The type for ray. This type is required only if
some dual functions are called.
*/
typedef unspecified_type Ray_2;

/*!
A function object to perform an incircle test for a point and three other points.
Provides the operator:

`Oriented_side operator()(Point p, Point q, Point r, Point s)`
which takes four points `p, q, r, s` as arguments and returns
`ON_POSITIVE_SIDE`, `ON_NEGATIVE_SIDE` or,
`ON_ORIENTED_BOUNDARY` according to the position of points `s`
with respect to the oriented circle through `p, q` and `r`.
*/
typedef unspecified_type Side_of_oriented_circle_2;

/*!
A function object to compare two distances for three points.
Provides the operator:

`Comparison_result operator()(Point_2 p, Point_2 q, Point_2 r)`
which returns `SMALLER`, `EQUAL` or `LARGER`
according to the distance between `p` and `q` being smaller, equal or larger
than the distance between `p` and `r`. This type is only require if
`nearest_vertex` queries are issued.
*/
typedef unspecified_type Compare_distance_2;

/*!
A function object to construct the circumcenter of three points.
Provides the operator:

`Point_2 operator()(Point_2 p, Point_2 q, Point_2 r)` which returns
the circumcenter of the three points `p, q` and `r`.
This type is required only if functions
relative to the dual Voronoi diagram are called.
*/
typedef unspecified_type Construct_circumcenter_2;

/*!
A function object to construct the bisector of two points.

Provides the operator:

`Line_2 operator()(Point_2 p, Point_2 q)` which constructs the
bisector line of points `p` and `q`.
This type is required only if functions
relative to the dual Voronoi diagram are called.
*/
typedef unspecified_type Construct_bisector_2;

/*!
A function object to build a ray from a point and a line.
Provides the operator:

`Ray_2 operator() ( Point_2 p, Line_2 l);`
*/
typedef unspecified_type Construct_ray_2;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
DelaunayTriangulationTraits_2();

/*!
copy constructor
*/
DelaunayTriangulationTraits_2(DelaunayTriangulationTraits_2
dtt);

/*!
Assignment operator.
*/
DelaunayTriangulationTraits_2
operator=(traits2);

/// @}

/// \name Access to Predicate and Constructor Objects
/// @{

/*!

*/
Side_of_oriented_circle_2
side_of_oriented_circle_2_object();

/// @}

/// \name
/// The following functions are required only if member functions of
/// the Delaunay triangulation relative to the dual Voronoi diagram
/// are called.
/// @{

/*!

*/
Compare_distance_2
compare_distance_2_object();

/*!

*/
Construct_circumcenter_2 construct_circumcenter_2_object();

/*!

*/
Construct_bisector_2 construct_bisector_2_object();


/*!

*/
Construct_ray_2 construct_ray_2_object();

/// @}

}; /* end DelaunayTriangulationTraits_2 */


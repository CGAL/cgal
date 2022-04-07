
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `ImplicitSurfaceTraits_3` describes the requirements of the traits class to
be plugged as `Traits` in `CGAL::Implicit_surface_3<Traits, Function>`.

When `make_surface_mesh` is called
with a surface of type `CGAL::Implicit_surface_3<Traits,Function>`,
the surface mesher traits generator generates automatically
a traits class that is a model of `SurfaceMeshTraits_3`.
Actually,
the concept `ImplicitSurfaceTraits_3` provides the types, predicates and constructors
that are passed
to the generated model of `SurfaceMeshTraits_3`.

\cgalHasModel Any \cgal Kernel.

\sa `CGAL::Implicit_surface_3<Traits, Function>`
\sa `CGAL::make_surface_mesh()`

*/

class ImplicitSurfaceTraits_3 {
public:

/// \name Types
/// @{

/*!
The numerical type. It must be model of
`FieldWithSqrt` and constructible from a `double`.
*/
typedef unspecified_type FT;

/*!
The point type. This point type must have a
constructor `Point_3(FT, FT, FT)`.
*/
typedef unspecified_type Point_3;

/*!
The line type.
*/
typedef unspecified_type Line_3;

/*!
The ray type.
*/
typedef unspecified_type Ray_3;

/*!
The segment type.
*/
typedef unspecified_type Segment_3;

/*!
The vector type.
*/
typedef unspecified_type Vector_3;

/*!
The sphere type.
*/
typedef unspecified_type Sphere_3;

/*!
A function object providing the operator

`FT operator()(Vector_3 v, Vector_3 w)` which returns the scalar
(inner) product of the two vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
A function object providing the operator

`FT operator()(Point_3, Point_3)` which returns the squared distance
between two points.
*/
typedef unspecified_type Compute_squared_distance_3;

/*!
A function object providing the operator

`FT operator()(const Sphere_3& s)` which returns the squared radius
of `s`.
*/
typedef unspecified_type Compute_squared_radius_3;

/*!
A function object providing the operator

`Point_3 operator()(const Sphere_3& s)` which computes the center of
the sphere `s`.
*/
typedef unspecified_type Construct_center_3;

/*!
A function object providing the operator

`Point_3 operator()(const Point_3& p, const Point_3& q)` which computes
the midpoint of the segment `pq`.
*/
typedef unspecified_type Construct_midpoint_3;

/*!
A function object providing the following operators:

`Point_3 operator()(const Line_3& l,int i)` which returns an
arbitrary point on `l`. It holds `point(i) == point(j)`, iff
`i==j`. Furthermore, is directed from `point(i)` to
`point(j)`, for all `i` \f$ <\f$ `j`.

`Point_3 operator()(const Ray_3& r,int i)` which returns a point on
`r`. `point(0)` is the source, `point(i)`, with
\f$ i>0\f$, is different from the source. \pre \f$ i \geq0\f$.

`Point_3 operator()(const Segment_3& s,int i)` which returns source
or target of `s`: `point(0)` returns the source of `s`,
`point(1)` returns the target of `s`. The parameter
`i` is taken modulo 2, which gives easy access to the other end
point.

*/
typedef unspecified_type Construct_point_on_3;

/*!
A function object providing the operator

`Segment_3 operator()(const Point_3 &p, const Point_3 &q)` which
returns a segment with source `p` and target `q`. It is directed from the
source towards the target.
*/
typedef unspecified_type Construct_segment_3;

/*!
A function object providing the operator

`Vector_3 operator()(const Vector_3 &v, const FT& scale)` which returns
the vector `v` scaled by a factor `scale`.
*/
typedef unspecified_type Construct_scaled_vector_3;

/*!
A function object providing the operator

`Point_3 operator()(const Point_3& p, const Vector_3& v)` which returns
the point obtained by translating `p` by the vector `v`.
*/
typedef unspecified_type Construct_translated_point_3;

/*!
A function object providing the operator

`Vector_3 operator()(const Point_3 &a, const Point_3 &b)` which returns
the vector `b-a`.
*/
typedef unspecified_type Construct_vector_3;

/*!
A function object providing the operator

`bool operator()(const Sphere_3&s, const Point_3&p)` which
returns true iff `p` lies on the bounded side of `s`.
*/
typedef unspecified_type Has_on_bounded_side_3;

/// @}

/// \name Operations
/// The following functions give access to the predicate and
/// construction objects:
/// @{

/*!

*/
Compute_scalar_product_3 compute_scalar_product_3_object();

/*!

*/
Compute_squared_distance_3 compute_squared_distance_3_object();

/*!

*/
Compute_squared_radius_3 compute_squared_radius_3_object();

/*!

*/
Construct_center_3 construct_center_3_object();

/*!

*/
Construct_midpoint_3 construct_midpoint_3_object();

/*!

*/
Construct_point_on_3 construct_point_on_3_object();

/*!

*/
Construct_scaled_vector_3 construct_scaled_vector_3_object();

/*!

*/
Construct_segment_3 construct_segment_3_object();

/*!

*/
Construct_translated_point_3 construct_translated_point_3_object();

/*!

*/
Construct_vector_3 construct_vector_3_object();

/*!

*/
Has_on_bounded_side_3 has_on_bounded_side_3_object();

/// @}

}; /* end ImplicitSurfaceTraits_3 */


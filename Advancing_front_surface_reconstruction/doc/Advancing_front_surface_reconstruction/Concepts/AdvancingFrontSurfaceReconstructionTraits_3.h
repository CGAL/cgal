
/*!
\ingroup PkgAdvancingFrontSurfaceReconstructionRef
\cgalConcept

The concept `AdvancingFrontSurfaceReconstructionTraits_3` describes the requirements
for the the geometric traits of the class `CGAL::Delaunay_triangulation_3`
used in the class `CGAL::Advancing_front_surface_reconstruction`.
It defines the geometric objects (points, segments...) forming the triangulation
together with a few geometric predicates and constructions on these objects.

\cgalRefines `DelaunayTriangulationTraits_3`

\cgalHasModel All models of `Kernel`.
*/
class AdvancingFrontSurfaceReconstructionTraits_3
{
public:

/// \name Types
/// @{

/*!
The coordinate type.
*/
typedef unspecified_type FT;

/*!
The vector type.
*/
typedef unspecified_type Vector_3;

/*!
The sphere type.
*/
typedef unspecified_type Sphere_3;

/*!
A constructor object that must provide the function operator

`Vector_3 operator()(Point_3 p, Point_3 q)`,

which constructs the vector `q-p`.
*/
typedef unspecified_type Construct_vector_3;

/*!
A constructor object that must provide the function operator

`Vector_3 operator()(Vector_3 v, Vector_3 w)`,

which returns the cross product of `v` and `w`.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/*!
A constructor object that must provide the function operator

`FT operator()(Vector_3 v, Vector_3 w)`,

which returns the scalar (inner) product of `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
A constructor object that must provide the function operator

`Sphere_3 operator()(Point_3 p, Point_3 q, Point_3 r)`,

which constructs a sphere initialized to the smallest sphere which passes
through the points `p`, `q`, and `r`.
*/
typedef unspecified_type Construct_sphere_3;

/*!
A constructor object that must provide the function operator

`Point_3 operator()(Sphere_3 s)`,

which returns the center of the sphere `s`.
*/
typedef unspecified_type Construct_center_3;

/*!
A constructor object that must provide the function operators

`FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`,

which returns the squared radius of the sphere passing through `p`, `q` and `r`,
and whose center is in the plane defined by these three points.

and

`FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`,

which returns the squared radius of the sphere passing through `p`, `q`, `r`, and `s`.

and

`FT operator()(Sphere_3 s)`,

which returns the squared radius of the sphere `s`.
*/
typedef unspecified_type Compute_squared_radius_3;

/*!
A constructor object that must provide the function operator

`FT operator()(Point_3 p, Point_3 q)`,

which returns the squared distance between the points `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_3;

/// @}

/// \name Operations
/// The following functions give access to the predicate and construction objects:
/// @{

/*!
gives access to the `Construct_vector_3` construction.
*/
Construct_vector_3 construct_vector_3_object();

/*!
gives access to the `Construct_cross_product_vector_3` construction.
*/
Construct_cross_product_vector_3 construct_cross_product_vector_3_object();

/*!
gives access to the `Compute_scalar_product_3` construction.
*/
Compute_scalar_product_3 compute_scalar_product_3_object();

/*!
gives access to the `Construct_sphere_3` construction.
*/
Construct_sphere_3 construct_sphere_3_object();

/*!
gives access to the `Construct_center_3` construction.
*/
Construct_center_3 construct_center_3_object();

/*!
gives access to the `Compute_squared_radius_3` construction.
*/
Compute_squared_radius_3 compute_squared_radius_3_object();

/*!
gives access to the `Compute_squared_distance_3` construction.
*/
Compute_squared_distance_3 compute_squared_distance_3_object();

/// @}

}; /* end AdvancingFrontSurfaceReconstructionTraits_3 */


/*!
\ingroup PkgMesh_3SecondaryConcepts
\cgalConcept

The concept `BisectionGeometricTraits_3` describes a geometric traits class 
that provides the basic types and operations 
to implement a model of `MeshDomain_3` 
based solely on intersection detections. 
Points in the non-empty intersections are herein computed 
by bisection. 

Such traits class is relevant when intersection detections 
can be performed efficiently. For instance, when bounding surfaces 
are implicitly described by a function (such as an isosurface of a 3D 
function from \f$ \mathbb{R}^3\f$ to \f$ \mathbb{R}\f$), the do-intersect predicate with a segment 
is computed by evaluations of the function values at both end points 
of the segment. 

\cgalHasModel Any \cgal Kernel.

\sa `ImplicitSurfaceTraits_3` 
\sa `IntersectionGeometricTraits_3` 
\sa `CGAL::Implicit_mesh_domain_3<Function,BGT>` 
\sa `CGAL::Labeled_image_mesh_domain_3<Image,BGT>` 

*/

class BisectionGeometricTraits_3 {
public:

/// \name Types 
/// @{

/*!
Numerical type. Must be a model of `::FieldNumberType` and
`::FieldWithSqrt`, and constructible from a `double`. 
*/ 
typedef unspecified_type FT; 

/*!
The point type. Must have a 
constructor `Point_3(FT, FT, FT)`. 
*/ 
typedef unspecified_type Point_3; 

/*!
Segment type. 
*/ 
typedef unspecified_type Segment_3; 

/*!
Ray type. 
*/ 
typedef unspecified_type Ray_3; 

/*!
Line type. 
*/ 
typedef unspecified_type Line_3; 

/*!
Vector type. 
*/ 
typedef unspecified_type Vector_3; 

/*!
Sphere type. 
*/ 
typedef unspecified_type Sphere_3; 

/*!
Model of `::Kernel::ComputeScalarProduct_3`.

That function object must provide the operator:
- `FT operator()(Vector_3 v, Vector_3 w)` which returns the scalar 
  (inner) product of the two vectors `v` and `w`. 
*/ 
typedef unspecified_type Compute_scalar_product_3; 

/*!
Model of `::Kernel::ComputeSquaredDistance_3`.

That function object must provide the operator:
- `FT operator()(Point_3, Point_3)` which returns the squared distance 
between two points. 
*/ 
typedef unspecified_type Compute_squared_distance_3; 

/*!
Model of `::Kernel::ComputeSquaredRadius_3`.

That function object must provide the operator:
- `FT operator()(Sphere_3 s)` which returns the squared radius 
of `s`. 
*/ 
typedef unspecified_type Compute_squared_radius_3; 

/*!
Model of `::Kernel::ConstructCenter_3`.

That function object must provide the operator:
- `Point_3 operator()(Sphere_3 s)` which returns the center of 
the sphere `s`. 
*/ 
typedef unspecified_type Construct_center_3; 

/*!
Model of `::Kernel::ConstructMidpoint_3`.

That function object must provide the operator:
- `Point_3 operator()(Point_3 p, Point_3 q)` which computes 
the midpoint of the segment `pq`. 
*/ 
typedef unspecified_type Construct_midpoint_3; 

/*!
Model of `::Kernel::ConstructPoint_3`.

That function object must provide the following operators:
- `Point_3 operator()(Line_3 l,int i)` which returns an 
arbitrary point on `l`. It holds `point(i) == point(j)`, iff 
`i==j`. Furthermore, is directed from `point(i)` to 
`point(j)`, for all `i` < `j`. 
- `Point_3 operator()(Ray_3 r,int i)` which returns a point on 
`r`. `point(0)` is the source, `point(i)`, with 
\f$ i>0\f$, is different from the source. \pre \f$ i \geq0\f$. 
- `Point_3 operator()(Segment_3 s,int i)` which returns either source 
or target of `s`: `point(0)` returns the source of `s`, 
`point(1)` returns the target of `s`. Parameter 
`i` is taken modulo 2, which gives easy access to the other end 
point. 

*/ 
typedef unspecified_type Construct_point_on_3; 

/*!
Model of `::Kernel::ConstructSegment_3`.

That function object must provide the operator:
- `Segment_3 operator()(Point_3 p, Point_3 q)` which 
returns a segment with source `p` and target `q`, directed from the 
source to the target. 
*/ 
typedef unspecified_type Construct_segment_3; 

/*!
Model of `::Kernel::ConstructScaledVector_3``.

That function object must provide the operator:
- `Vector_3 operator()(Vector_3 v, FT scale)` which returns 
the vector `v` scaled by a factor `scale`. 
*/ 
typedef unspecified_type Construct_scaled_vector_3; 

/*!
Model of `::Kernel::ConstructTranslatedPoint_3`.

That function object must provide the operator:
- `Point_3 operator()(Point_3 p, Vector_3 v)` which returns 
the point obtained by translating `p` by the vector `v`. 
*/ 
typedef unspecified_type Construct_translated_point_3; 

/*!
Model of `::Kernel::ConstructVertex_3`.

That function object must provide the operator:
- `Vector_3 operator()(Point_3 a, Point_3 b)` which returns 
the vector `b-a`. 
*/ 
typedef unspecified_type Construct_vector_3; 

/*!
Model of `::Kernel::HasOnBoundedSide_3`.

That function object must provide the operator:
- `bool operator()(Sphere_3 s, Point_3 p)` which 
returns true iff `p` lies on the bounded side of `s. 
*/ 
typedef unspecified_type Has_on_bounded_side_3; 

/// @} 

/// \name Operations 
/// The following functions give access to the predicate and construction objects: 
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

}; /* end BisectionGeometricTraits_3 */

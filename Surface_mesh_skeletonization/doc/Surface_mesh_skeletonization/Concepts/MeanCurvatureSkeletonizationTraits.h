
/*!
 * \ingroup PkgSurfaceMeshSkeletonizationConcepts
 * \cgalConcept
 *
 * Traits class concept defining the requirements of the class `CGAL::Mean_curvature_flow_skeletonization`.
 *
 * \cgalHasModel Any \cgal `Kernel` with `double` as `%Kernel::%FT`
 *
 *
 */
class MeanCurvatureSkeletonizationTraits {
public:

/// \name Types
/// @{

/// The point type.
typedef unspecified_type Point_3;


/// The vector type.
typedef unspecified_type Vector_3;

/// The number type. It must be constructible from and convertible to `double`.
typedef unspecified_type FT;

/*!
 * Function object type that provides
 * `Point_3 operator()(FT x, FT y, FT z) const`
 * returning the point with `x`, `y` and `z` as Cartesian coordinates.
 */
typedef unspecified_type Construct_point_3;

/*!
 * Function object type that provides `Point_3 operator()(Point_3 p1, Point_3 p2) const`
 * returning the vector `p1p2`, and `Point_3 operator()(NULL_VECTOR) const` returning the null vector.
 */
typedef unspecified_type Construct_vector_3;

/*!
 * Function object type that provides
 * `Vector_3 operator()(Vector_3 v, FT t) const`
 * returning the vector `t * v`.
 */
typedef unspecified_type Construct_scaled_vector_3;

/*!
 * Function object type that provides
 * `Vector_3 operator()(Vector_3 v, FT t) const`
 * returning the vector `v / t`.
 */
typedef unspecified_type Construct_divided_vector_3;

/*!
 * Function object type that provides
 * `Vector_3 operator()(Vector_3 v1, Vector_3 v2) const`
 * returning the cross-product vector of `v1` and `v2`.
 */
typedef unspecified_type Construct_cross_product_vector_3;

/*!
 * Function object type that provides
 * `Vector_3 operator()(Vector_3 v1, Vector_3 v2) const`
 * returning the `v1+v2`.
 */
typedef unspecified_type Construct_sum_of_vectors_3;

/*!
 * Function object type that provides
 * `Point_3 operator()(Point_3 p1, Point_3 p2) const`
 * returning the midpoint of `p1` and `p2`.
 */
typedef unspecified_type Construct_midpoint_3;

/*!
 * Function object type that provides
 * `FT operator()(Point_3 p, Point_3 q) const`
 * returning the squared distance between `p` and `q`.
 */
typedef unspecified_type Compute_squared_distance_3;

/*!
 * Function object type that provides
 * `FT operator()(Vector_3 v) const`
 * returning the squared length of `v`.
 */
typedef unspecified_type Compute_squared_length_3;

/*!
 * Function object type that provides
 * `FT operator()(Point_3 p1,Point_3 p2,Point_3 p3) const`
 * returning the area of the triangle defined by `p1`, `p2` and `p3`.
 */
typedef unspecified_type Compute_area_3;

/*!
 * Function object type that provides
 * `FT operator()(Vector_3 v1, Vector_3 v2) const`
 * returning the scalar product of `v1` and `v2`.
 */
typedef unspecified_type Compute_scalar_product_3;

/*!
 * Function object type that provides
 * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v) const`
 * returning the `x` coordinate of a point and a vector respectively.
 */
typedef unspecified_type Compute_x_3;

/*!
 * Function object type that provides
 * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v) const`
 * returning the `y` coordinate of a point and a vector respectively.
 */
typedef unspecified_type Compute_y_3;

/*!
 * Function object type that provides
 * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v) const`
 * returning the `z` coordinate of a point and a vector respectively.
 */
typedef unspecified_type Compute_z_3;

/// @}

/// \name Access to Function Objects
/// @{

/// Function object creator
Construct_point_3
construct_point_3_object();

/// Function object creator
Construct_vector_3
construct_vector_3_object();

/// Function object creator
Construct_scaled_vector_3
construct_scaled_vector_3_object();

/// Function object creator
Construct_divided_vector_3
construct_divided_vector_3_object();

/// Function object creator
Construct_cross_product_vector_3
construct_cross_product_vector_3_object();

/// Function object creator
Construct_sum_of_vectors_3
construct_sum_of_vectors_3_object();

/// Function object creator
Construct_midpoint_3
construct_midpoint_3_object();

/// Function object creator
Compute_squared_distance_3
compute_squared_distance_3_object();

/// Function object creator
Compute_squared_length_3
compute_squared_length_3_object();

/// Function object creator
Compute_area_3
compute_area_3_object();

/// Function object creator
Compute_scalar_product_3
compute_scalar_product_3_object();

/// Function object creator
Compute_x_3
compute_x_3_object();

/// Function object creator
Compute_y_3
compute_y_3_object();

/// Function object creator
Compute_z_3
compute_z_3_object();

/// @}

}; /* end DelaunayTriangulationTraits_2 */

